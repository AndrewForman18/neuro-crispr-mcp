"""NeuroPlex Base Ingestor — Abstract framework for source-specific ingestion.

Each external source implements a subclass of BaseIngestor that handles:
  1. API authentication and rate limiting
  2. Data fetching (paginated or bulk)
  3. Schema normalization to the common base + VARIANT payload
  4. Delta table upsert (MERGE for incremental, overwrite for bulk)

Usage from a Databricks notebook:
    from ingestion.source_registry import SOURCE_MAP
    from ingestion.ingestors.opentargets import OpenTargetsIngestor

    ingestor = OpenTargetsIngestor(SOURCE_MAP["opentargets"])
    ingestor.run(gene="HCRT")
"""

from __future__ import annotations

import json
import logging
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Generator, Optional

import requests

from .source_registry import SourceConfig, IngestMethod, TARGET_CATALOG, TARGET_SCHEMA

logger = logging.getLogger("neuroplex.ingestor")


# ═══════════════════════════════════════════════════════════════════════════════
# Rate Limiter
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class RateLimiter:
    """Simple token-bucket rate limiter."""
    requests_per_minute: int
    _timestamps: list[float] = field(default_factory=list, repr=False)

    def wait(self):
        """Block until a request slot is available."""
        now = time.time()
        window_start = now - 60.0
        # Purge old timestamps
        self._timestamps = [t for t in self._timestamps if t > window_start]
        if len(self._timestamps) >= self.requests_per_minute:
            sleep_time = self._timestamps[0] - window_start + 0.1
            logger.debug(f"Rate limit: sleeping {sleep_time:.1f}s")
            time.sleep(sleep_time)
        self._timestamps.append(time.time())


# ═══════════════════════════════════════════════════════════════════════════════
# Normalized Record
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class NormalizedRecord:
    """Standardized record matching the target table schema."""
    record_id: str
    source_key: str
    payload: dict[str, Any]
    gene_symbol: Optional[str] = None
    disease: Optional[str] = None
    drug: Optional[str] = None
    title: Optional[str] = None
    summary: Optional[str] = None
    source_updated_at: Optional[datetime] = None

    def to_row(self) -> dict[str, Any]:
        """Convert to a dict suitable for Spark DataFrame creation."""
        return {
            "record_id": self.record_id,
            "source_key": self.source_key,
            "gene_symbol": self.gene_symbol,
            "disease": self.disease,
            "drug": self.drug,
            "title": self.title,
            "summary": self.summary,
            "payload": json.dumps(self.payload),
            "ingested_at": datetime.now(timezone.utc).isoformat(),
            "source_updated_at": self.source_updated_at.isoformat() if self.source_updated_at else None,
        }


# ═══════════════════════════════════════════════════════════════════════════════
# Base Ingestor (Abstract)
# ═══════════════════════════════════════════════════════════════════════════════

class BaseIngestor(ABC):
    """Abstract base class for all NeuroPlex data source ingestors."""

    def __init__(self, config: SourceConfig):
        self.config = config
        self.rate_limiter = RateLimiter(config.rate_limit_rpm)
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "NeuroPlex/1.0 (Eisai Discovery)"})
        self._setup_auth()

    def _setup_auth(self):
        """Configure authentication if required."""
        if self.config.auth_required and self.config.auth_secret_key:
            try:
                from databricks.sdk import WorkspaceClient
                w = WorkspaceClient()
                # Assume secrets are in 'neuroplex' scope
                secret = w.secrets.get_secret("neuroplex", self.config.auth_secret_key)
                self._apply_auth(secret.value)
            except Exception as e:
                logger.warning(f"Auth setup failed for {self.config.name}: {e}")

    def _apply_auth(self, secret_value: str):
        """Apply auth credentials to session. Override for custom auth patterns."""
        # Default: add as API key query parameter
        # Subclasses override for Bearer tokens, custom headers, etc.
        self._api_key = secret_value

    def _request(self, method: str, url: str, **kwargs) -> requests.Response:
        """Rate-limited HTTP request with retry logic."""
        self.rate_limiter.wait()
        retries = 3
        for attempt in range(retries):
            try:
                resp = self.session.request(method, url, timeout=30, **kwargs)
                if resp.status_code == 429:  # Rate limited
                    retry_after = int(resp.headers.get("Retry-After", 60))
                    logger.warning(f"429 from {self.config.name}, sleeping {retry_after}s")
                    time.sleep(retry_after)
                    continue
                resp.raise_for_status()
                return resp
            except requests.exceptions.RequestException as e:
                if attempt < retries - 1:
                    time.sleep(2 ** attempt)
                    continue
                raise
        raise RuntimeError(f"Max retries exceeded for {url}")

    def get(self, url: str, **kwargs) -> requests.Response:
        return self._request("GET", url, **kwargs)

    def post(self, url: str, **kwargs) -> requests.Response:
        return self._request("POST", url, **kwargs)

    # ───── Abstract methods (implement per source) ─────

    @abstractmethod
    def fetch(self, **query_params) -> Generator[dict[str, Any], None, None]:
        """Fetch raw records from the source API.

        Yields raw API response records (dicts) one at a time.
        Handles pagination internally.
        """
        ...

    @abstractmethod
    def normalize(self, raw_record: dict[str, Any]) -> NormalizedRecord:
        """Transform a raw API record into a NormalizedRecord.

        Maps source-specific fields to the common schema.
        """
        ...

    # ───── Orchestration ─────

    def run(
        self,
        spark=None,
        mode: str = "merge",
        batch_size: int = 500,
        **query_params,
    ) -> dict[str, Any]:
        """Full ingestion pipeline: fetch → normalize → write to Delta.

        Args:
            spark: SparkSession (auto-detected if None)
            mode: 'merge' for incremental upsert, 'overwrite' for full refresh
            batch_size: Records per write batch
            **query_params: Source-specific query parameters

        Returns:
            Summary dict with record counts and timing.
        """
        if spark is None:
            from pyspark.sql import SparkSession
            spark = SparkSession.getActiveSession()

        start = time.time()
        fqn = f"{self.config.target_table}"
        full_fqn = f"{TARGET_CATALOG}.{TARGET_SCHEMA}.{fqn}"

        logger.info(f"Ingesting {self.config.name} → {full_fqn} (mode={mode})")

        records = []
        errors = 0
        for raw in self.fetch(**query_params):
            try:
                normalized = self.normalize(raw)
                records.append(normalized.to_row())
            except Exception as e:
                errors += 1
                logger.warning(f"Normalize error: {e}")
                continue

            # Batch write
            if len(records) >= batch_size:
                self._write_batch(spark, records, full_fqn, mode)
                records = []
                mode = "merge"  # After first batch, always merge

        # Final batch
        if records:
            self._write_batch(spark, records, full_fqn, mode)

        elapsed = time.time() - start
        summary = {
            "source": self.config.name,
            "records_written": len(records),
            "errors": errors,
            "elapsed_seconds": round(elapsed, 1),
            "target_table": full_fqn,
        }
        logger.info(f"Ingestion complete: {summary}")
        return summary

    def _write_batch(self, spark, rows: list[dict], table_fqn: str, mode: str):
        """Write a batch of normalized records to Delta."""
        from pyspark.sql import functions as F

        df = spark.createDataFrame(rows)

        if mode == "overwrite":
            df.write.mode("overwrite").saveAsTable(table_fqn)
        else:
            # MERGE on record_id + source_key
            df.createOrReplaceTempView("_neuroplex_batch")
            spark.sql(f"""
                MERGE INTO {table_fqn} AS target
                USING _neuroplex_batch AS source
                ON target.record_id = source.record_id
                   AND target.source_key = source.source_key
                WHEN MATCHED THEN UPDATE SET *
                WHEN NOT MATCHED THEN INSERT *
            """)
            spark.catalog.dropTempView("_neuroplex_batch")


# ═══════════════════════════════════════════════════════════════════════════════
# Query Executor (for MCP tool calls at runtime)
# ═══════════════════════════════════════════════════════════════════════════════

def query_source(source_key: str, query: str, limit: int = 20, **filters) -> str:
    """Generic query function for MCP tool calls.

    Queries the pre-ingested Delta table for a source.
    Used by the dynamic tool executor in app.py.

    Args:
        source_key: Registry key (e.g., 'opentargets', 'gnomad')
        query: Search term (gene, drug, disease, free text)
        limit: Max results
        **filters: Additional field-specific filters

    Returns:
        JSON string of results (for LLM consumption)
    """
    from .source_registry import SOURCE_MAP, TARGET_CATALOG, TARGET_SCHEMA

    source = SOURCE_MAP.get(source_key)
    if not source:
        return json.dumps({"error": f"Unknown source: {source_key}"})

    table = f"{TARGET_CATALOG}.{TARGET_SCHEMA}.{source.target_table}"

    # Build WHERE clause
    conditions = []
    params = []

    # Primary search: match across gene_symbol, disease, drug, title, summary
    search_clause = " OR ".join([
        f"LOWER(gene_symbol) LIKE LOWER('%{query}%')",
        f"LOWER(disease) LIKE LOWER('%{query}%')",
        f"LOWER(drug) LIKE LOWER('%{query}%')",
        f"LOWER(title) LIKE LOWER('%{query}%')",
        f"LOWER(summary) LIKE LOWER('%{query}%')",
    ])
    conditions.append(f"({search_clause})")

    # Additional filters
    for field_name, value in filters.items():
        if value and field_name in source.searchable_fields:
            conditions.append(f"LOWER({field_name}) LIKE LOWER('%{value}%')")

    where = " AND ".join(conditions) if conditions else "1=1"
    sql = f"""
        SELECT record_id, gene_symbol, disease, drug, title, summary,
               payload, source_updated_at
        FROM {table}
        WHERE {where}
        ORDER BY source_updated_at DESC NULLS LAST
        LIMIT {limit}
    """

    # Execute via the existing SQL connection
    try:
        # Import the shared executor from the MCP server
        import sys, os
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        from neuro_mcp_server.server import _execute_query
        results = _execute_query(sql, timeout=30)
        return json.dumps({"source": source.name, "query": query, "count": len(results), "results": results})
    except Exception as e:
        return json.dumps({"error": str(e), "source": source.name, "query": query})

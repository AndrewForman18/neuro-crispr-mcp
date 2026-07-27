"""NCBI GTR Ingestor — Genetic Testing Registry for gene-based diagnostic tests.

Uses NCBI E-utilities API:
  - esearch.fcgi  — search GTR by gene name, returns GTR test IDs
  - esummary.fcgi — batch-fetch full test metadata (lab, conditions, methodology)

Field mapping confirmed from live GTR esummary response (July 2025):
  testname, accession, testtype, analytes, conditionlist,
  offerer, offererlocation, country, method, testpurpose

API: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
GTR: https://www.ncbi.nlm.nih.gov/gtr/
"""
from __future__ import annotations

import logging
from typing import Any, Generator

from ..base_ingestor import BaseIngestor, NormalizedRecord
from ..source_registry import SourceConfig

logger = logging.getLogger("neuroplex.ncbi_gtr")

_ESEARCH  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
_BATCH    = 200   # max IDs per esummary call


class NCBIGTRIngestor(BaseIngestor):
    """Ingest genetic test records from NCBI Genetic Testing Registry (GTR)."""

    def _apply_auth(self, secret_value: str):
        """Store NCBI API key for higher rate limits (10 req/s vs 3/s)."""
        self._api_key = secret_value

    # ── Helpers ────────────────────────────────────────────────────────────────

    def _params(self, extra: dict | None = None) -> dict:
        """Build base params, injecting NCBI API key when available."""
        p: dict = {}
        key = getattr(self, "_api_key", None)
        if key:
            p["api_key"] = key
        if extra:
            p.update(extra)
        return p

    def _search(self, gene: str, retmax: int = 100) -> list[str]:
        """Return GTR test IDs for tests that target *gene* via [Gene_Name] field."""
        resp = self.get(_ESEARCH, params=self._params({
            "db":      "gtr",
            "term":    f'"{gene}"[Gene_Name]',
            "retmax":  retmax,
            "retmode": "json",
        }))
        return resp.json().get("esearchresult", {}).get("idlist", [])

    def _summaries(self, ids: list[str]) -> list[dict]:
        """Batch-fetch esummary records; returns a flat list of summary dicts."""
        results: list[dict] = []
        for i in range(0, len(ids), _BATCH):
            chunk = ids[i : i + _BATCH]
            try:
                resp = self.get(_ESUMMARY, params=self._params({
                    "db":      "gtr",
                    "id":      ",".join(chunk),
                    "retmode": "json",
                }))
                data = resp.json().get("result", {})
                uids = data.get("uids", [])
                results.extend(data[uid] for uid in uids if uid in data)
            except Exception as e:
                logger.warning(f"GTR esummary batch failed for {len(chunk)} IDs: {e}")
        return results

    # ── BaseIngestor interface ──────────────────────────────────────────────────

    def fetch(self, gene: str = "", limit: int = 100, **kwargs) -> Generator[dict[str, Any], None, None]:
        """Fetch GTR test records for a gene.

        Args:
            gene:  HUGO gene symbol (e.g. HCRT, PSEN1, APP)
            limit: Maximum number of tests to return
        Yields:
            Raw esummary dicts, each augmented with _query_gene key.
        """
        if not gene:
            return
        ids = self._search(gene, retmax=limit)
        if not ids:
            logger.info(f"GTR: no tests found for {gene}")
            return
        logger.info(f"GTR: {len(ids)} test IDs for {gene}, fetching summaries…")
        for summary in self._summaries(ids):
            summary["_query_gene"] = gene.upper()
            yield summary

    def normalize(self, raw: dict[str, Any]) -> NormalizedRecord:
        """Normalise an NCBI GTR esummary record into the NeuroPlex common schema.

        Field names from live GTR esummary API (July 2025):
          testname      → test display name
          accession     → GTR accession string (e.g. GTR000568998)
          testtype      → e.g. 'Clinical'
          analytes      → [{analytetype, name, geneid, location}, ...]
          conditionlist → [{name, acronym, cui}, ...]
          offerer       → lab/organisation name (plain string)
          offererlocation → {city, state, postcode, country, ...}
          country       → top-level country string
          method        → [{name, categoriesstring, categorylist}, ...]
          testpurpose   → ['Screening', ...]
        """
        uid        = str(raw.get("uid", ""))
        query_gene = raw.get("_query_gene", "")

        # ── Test identity ──
        test_name  = (raw.get("testname") or "").strip()
        accession  = (raw.get("accession") or f"GTR{uid}").strip()
        test_type  = (raw.get("testtype") or "").strip()

        # ── Genes tested via analytes list ──
        analytes     = raw.get("analytes") or []
        genes_tested = [
            a["name"] for a in analytes
            if isinstance(a, dict) and a.get("analytetype") == "Gene" and a.get("name")
        ]
        primary_gene = genes_tested[0] if genes_tested else query_gene

        # ── Conditions / diseases ──
        conditions = [
            c["name"] for c in (raw.get("conditionlist") or [])
            if isinstance(c, dict) and c.get("name")
        ]

        # ── Lab / organisation ──
        lab_name    = (raw.get("offerer") or "").strip()
        loc         = raw.get("offererlocation") or {}
        lab_city    = (loc.get("city")  or "").strip()
        lab_state   = (loc.get("state") or "").strip()
        lab_country = (raw.get("country") or loc.get("country") or "").strip()
        location_parts = [p for p in [lab_city, lab_state, lab_country] if p]
        lab_location   = ", ".join(location_parts)

        # ── Methodology ──
        method_list  = raw.get("method") or []
        if isinstance(method_list, list):
            method = "; ".join(m["name"] for m in method_list if isinstance(m, dict) and m.get("name"))
        else:
            method = str(method_list)
        test_purpose = "; ".join(raw.get("testpurpose") or [])
        order_url    = raw.get("orderurl") or ""

        return NormalizedRecord(
            record_id   = f"gtr_{uid}",
            source_key  = "ncbi_gtr",
            gene_symbol = primary_gene or query_gene,
            disease     = conditions[0] if conditions else None,
            title       = f"{test_name} [{accession}]",
            summary     = (
                f"Test: {test_name or 'N/A'} | "
                f"Gene(s): {', '.join(genes_tested) or query_gene} | "
                f"Conditions: {', '.join(conditions[:3]) or 'Not specified'} | "
                f"Method: {method or 'Not specified'} | "
                f"Type: {test_type or 'Not specified'} | "
                f"Lab: {lab_name or 'Not specified'}"
                + (f" ({lab_location})" if lab_location else "")
            ),
            payload = {
                "gtr_id":             uid,
                "accession":          accession,
                "test_name":          test_name,
                "test_type":          test_type,
                "test_purpose":       test_purpose,
                "genes_tested":       genes_tested,
                "conditions":         conditions,
                "method":             method,
                "lab_name":           lab_name,
                "lab_country":        lab_country,
                "lab_state":          lab_state,
                "lab_city":           lab_city,
                "order_url":          order_url,
                "analytical_validity": raw.get("analyticalvalidity") or "",
                "clinical_utility":    raw.get("clinicalutility")    or "",
                "raw_summary":        {k: v for k, v in raw.items() if k != "_query_gene"},
            },
        )

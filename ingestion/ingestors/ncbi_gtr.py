"""NCBI GTR Ingestor — Genetic Testing Registry.

API: NCBI E-utilities (esearch + efetch for GTR database)
"""
from __future__ import annotations
import json
from typing import Any, Generator
from ..base_ingestor import BaseIngestor, NormalizedRecord
from ..source_registry import SourceConfig


class NcbiGtrIngestor(BaseIngestor):
    BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def _apply_auth(self, secret_value: str):
        self._api_key = secret_value

    def fetch(self, gene: str = "", limit: int = 20, **kwargs) -> Generator[dict[str, Any], None, None]:
        if not gene:
            return
        params = {"db": "gtr", "term": f"{gene}[gene]", "retmax": limit, "retmode": "json"}
        if hasattr(self, "_api_key") and self._api_key:
            params["api_key"] = self._api_key
        resp = self.get(f"{self.BASE}/esearch.fcgi", params=params)
        data = resp.json()
        ids = data.get("esearchresult", {}).get("idlist", [])
        for gtr_id in ids[:limit]:
            yield {"record_type": "test", "gtr_id": gtr_id, "gene_symbol": gene.upper()}

    def normalize(self, raw: dict[str, Any]) -> NormalizedRecord:
        symbol = raw["gene_symbol"]
        gtr_id = raw["gtr_id"]
        return NormalizedRecord(
            record_id=f"gtr_{gtr_id}",
            source_key="ncbi_gtr",
            gene_symbol=symbol,
            title=f"{symbol} — GTR Test ID {gtr_id}",
            summary=f"NCBI Genetic Testing Registry entry for {symbol}. GTR ID: {gtr_id}",
            payload={"gtr_id": gtr_id, "gene": symbol},
        )

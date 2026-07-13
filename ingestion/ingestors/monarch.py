"""Monarch Initiative Ingestor — Gene-disease associations and phenotype matching.

API: REST v3 at https://api.monarchinitiative.org/v3/api/
"""
from __future__ import annotations
import json
from typing import Any, Generator
from ..base_ingestor import BaseIngestor, NormalizedRecord
from ..source_registry import SourceConfig


class MonarchIngestor(BaseIngestor):
    BASE = "https://api.monarchinitiative.org/v3/api"

    def _apply_auth(self, secret_value: str):
        pass

    def fetch(self, gene: str = "", limit: int = 20, **kwargs) -> Generator[dict[str, Any], None, None]:
        if not gene:
            return
        # Search for gene entity
        resp = self.get(f"{self.BASE}/search", params={"q": gene, "category": "biolink:Gene", "limit": 5})
        results = resp.json().get("items", [])
        gene_id = None
        for item in results:
            if gene.upper() in (item.get("symbol") or "").upper():
                gene_id = item["id"]
                break
        if not gene_id and results:
            gene_id = results[0]["id"]
        if not gene_id:
            return

        # Get disease associations
        resp = self.get(f"{self.BASE}/entity/{gene_id}/associations",
                       params={"category": "biolink:GeneToDiseaseAssociation", "limit": limit})
        assocs = resp.json().get("items", [])
        for assoc in assocs:
            yield {"record_type": "disease_association", "assoc": assoc, "gene_symbol": gene.upper(), "gene_id": gene_id}

        # Get phenotype associations
        resp = self.get(f"{self.BASE}/entity/{gene_id}/associations",
                       params={"category": "biolink:GeneToPhenotypicFeatureAssociation", "limit": limit})
        phenos = resp.json().get("items", [])
        for pheno in phenos:
            yield {"record_type": "phenotype", "assoc": pheno, "gene_symbol": gene.upper(), "gene_id": gene_id}

    def normalize(self, raw: dict[str, Any]) -> NormalizedRecord:
        symbol = raw["gene_symbol"]
        assoc = raw["assoc"]
        obj = assoc.get("object", {})
        obj_name = obj.get("label", obj.get("id", "unknown"))

        if raw["record_type"] == "disease_association":
            return NormalizedRecord(
                record_id=f"monarch_{raw['gene_id']}_{obj.get('id', '')}",
                source_key="monarch",
                gene_symbol=symbol,
                disease=obj_name,
                title=f"{symbol} ↔ {obj_name}",
                summary=f"Disease association. Source: {', '.join(s.get('label', '') for s in assoc.get('publications', [])[:3])}. "
                        f"Evidence: {assoc.get('evidence_count', 'N/A')}",
                payload=assoc,
            )
        else:
            return NormalizedRecord(
                record_id=f"monarch_{raw['gene_id']}_pheno_{obj.get('id', '')}",
                source_key="monarch",
                gene_symbol=symbol,
                title=f"{symbol} → {obj_name} (phenotype)",
                summary=f"Phenotype association: {obj_name}. ID: {obj.get('id', '')}",
                payload=assoc,
            )

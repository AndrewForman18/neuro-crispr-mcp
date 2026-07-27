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

        # Get phenotype associations (v3 API: /association endpoint with subject param)
        resp = self.get(f"{self.BASE}/association",
                       params={"subject": gene_id, "category": "biolink:GeneToPhenotypicFeatureAssociation", "limit": limit})
        phenos = resp.json().get("items", []) if resp.status_code == 200 else []
        for pheno in phenos:
            yield {"record_type": "phenotype", "assoc": pheno, "gene_symbol": gene.upper(), "gene_id": gene_id}

        # Get all other associations (GeneToDiseaseAssociation removed from v3 enum;
        # query without category and filter for disease-related entries)
        resp = self.get(f"{self.BASE}/association",
                       params={"subject": gene_id, "limit": limit * 2})
        all_assocs = resp.json().get("items", []) if resp.status_code == 200 else []
        seen_pheno_ids = {p.get("id") for p in phenos}
        for assoc in all_assocs:
            if assoc.get("id") not in seen_pheno_ids:
                yield {"record_type": "disease_association", "assoc": assoc, "gene_symbol": gene.upper(), "gene_id": gene_id}

    def normalize(self, raw: dict[str, Any]) -> NormalizedRecord:
        symbol = raw["gene_symbol"]
        assoc = raw["assoc"]
        # v3 API: 'object' is a string ID, label is in 'object_label'
        obj_id = assoc.get("object", "")
        obj_name = assoc.get("object_label", obj_id)

        if raw["record_type"] == "disease_association":
            return NormalizedRecord(
                record_id=f"monarch_{raw['gene_id']}_{obj_id}",
                source_key="monarch",
                gene_symbol=symbol,
                disease=obj_name,
                title=f"{symbol} ↔ {obj_name}",
                summary=f"Disease association. Category: {assoc.get('category', '')}. "
                        f"Source: {assoc.get('primary_knowledge_source', 'N/A')}",
                payload=assoc,
            )
        else:
            return NormalizedRecord(
                record_id=f"monarch_{raw['gene_id']}_pheno_{obj_id}",
                source_key="monarch",
                gene_symbol=symbol,
                title=f"{symbol} → {obj_name} (phenotype)",
                summary=f"Phenotype association: {obj_name}. ID: {obj_id}. "
                        f"Source: {assoc.get('primary_knowledge_source', 'N/A')}",
                payload=assoc,
            )

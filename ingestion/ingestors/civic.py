"""CIViC Ingestor — Clinical variant evidence and therapy associations.

CIViC (Clinical Interpretation of Variants in Cancer):
- Curated clinical evidence for gene variants
- Therapy associations with evidence levels (A-E)
- Disease-variant-drug relationships

API: GraphQL at https://civicdb.org/api/graphql
"""
from __future__ import annotations
import json
from typing import Any, Generator
from ..base_ingestor import BaseIngestor, NormalizedRecord
from ..source_registry import SourceConfig

_GENE_QUERY = """
query GeneVariants($name: String!) {
  gene(name: $name) {
    id
    name
    description
    variants(first: 50) {
      nodes {
        id
        name
        singleVariantMolecularProfile {
          molecularProfileScore
        }
        variantTypes {
          name
        }
        evidenceItems(first: 20) {
          nodes {
            id
            status
            evidenceLevel
            evidenceType
            evidenceDirection
            significance
            disease { name doid }
            therapies { name ncitId }
            description
          }
        }
      }
    }
  }
}
"""


class CivicIngestor(BaseIngestor):
    CIVIC_API = "https://civicdb.org/api/graphql"

    def _apply_auth(self, secret_value: str):
        pass  # Open API

    def fetch(self, gene: str = "", limit: int = 50, **kwargs) -> Generator[dict[str, Any], None, None]:
        if not gene:
            return
        resp = self.post(self.CIVIC_API, json={"query": _GENE_QUERY, "variables": {"name": gene.upper()}})
        data = resp.json().get("data", {}).get("gene")
        if not data:
            return
        for variant in (data.get("variants") or {}).get("nodes", []):
            yield {"record_type": "variant", "gene": data, "variant": variant}
            for ev in (variant.get("evidenceItems") or {}).get("nodes", []):
                yield {"record_type": "evidence", "gene": data, "variant": variant, "evidence": ev}

    def normalize(self, raw: dict[str, Any]) -> NormalizedRecord:
        gene = raw["gene"]
        symbol = gene["name"]
        if raw["record_type"] == "variant":
            v = raw["variant"]
            return NormalizedRecord(
                record_id=f"civic_variant_{v['id']}",
                source_key="civic",
                gene_symbol=symbol,
                title=f"{symbol} {v['name']}",
                summary=f"Variant types: {', '.join(t['name'] for t in v.get('variantTypes', []))}. "
                        f"MP score: {v.get('singleVariantMolecularProfile', {}).get('molecularProfileScore', 'N/A')}",
                payload=v,
            )
        else:  # evidence
            ev = raw["evidence"]
            v = raw["variant"]
            disease = ev.get("disease", {}).get("name", "")
            therapies = ", ".join(t["name"] for t in ev.get("therapies", []))
            return NormalizedRecord(
                record_id=f"civic_evidence_{ev['id']}",
                source_key="civic",
                gene_symbol=symbol,
                disease=disease,
                drug=therapies or None,
                title=f"{symbol} {v['name']} — Level {ev.get('evidenceLevel', '?')} ({ev.get('evidenceType', '')})",
                summary=f"{ev.get('significance', '')}. Disease: {disease}. Therapies: {therapies or 'none'}. {(ev.get('description') or '')[:200]}",
                payload=ev,
            )

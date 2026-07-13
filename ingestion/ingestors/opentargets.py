"""OpenTargets Ingestor — Target-disease associations and druggability.

OpenTargets Platform GraphQL API:
- Target (gene) → disease associations with evidence scores
- Tractability (druggability) data
- Known drugs in clinical/approved stages
- Disease-gene evidence from multiple sources

API: https://api.platform.opentargets.org/api/v4/graphql
Docs: https://platform-docs.opentargets.org/
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from typing import Any, Generator

from ..base_ingestor import BaseIngestor, NormalizedRecord
from ..source_registry import SourceConfig


# GraphQL queries for OpenTargets
_TARGET_ASSOCIATIONS_QUERY = """
query TargetAssociations($ensemblId: String!, $size: Int!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    approvedName
    biotype
    tractability {
      label
      modality
      value
    }
    associatedDiseases(page: {size: $size, index: 0}) {
      count
      rows {
        disease {
          id
          name
          therapeuticAreas {
            id
            name
          }
        }
        score
        datasourceScores {
          id
          score
        }
      }
    }
    knownDrugs(size: $size) {
      count
      rows {
        drug {
          id
          name
          drugType
          maximumClinicalTrialPhase
          mechanismsOfAction {
            rows {
              mechanismOfAction
              actionType
            }
          }
        }
        disease {
          id
          name
        }
        phase
        status
        urls {
          niceName
          url
        }
      }
    }
  }
}
"""

_SEARCH_QUERY = """
query SearchTarget($queryString: String!, $size: Int!) {
  search(queryString: $queryString, entityNames: ["target"], page: {size: $size, index: 0}) {
    total
    hits {
      id
      entity
      name
      description
      ... on Target {
        approvedSymbol
      }
    }
  }
}
"""


class OpenTargetsIngestor(BaseIngestor):
    """Ingest target-disease association and druggability data from OpenTargets."""

    def _apply_auth(self, secret_value: str):
        """OpenTargets API is open — no auth needed."""
        pass

    def _graphql(self, query: str, variables: dict) -> dict:
        """Execute a GraphQL query against OpenTargets."""
        resp = self.post(
            self.config.base_url,
            json={"query": query, "variables": variables},
        )
        data = resp.json()
        if "errors" in data:
            raise RuntimeError(f"GraphQL errors: {data['errors']}")
        return data.get("data", {})

    def _resolve_ensembl_id(self, gene_symbol: str) -> str | None:
        """Resolve a gene symbol to an Ensembl ID via OpenTargets search."""
        data = self._graphql(_SEARCH_QUERY, {"queryString": gene_symbol, "size": 5})
        hits = data.get("search", {}).get("hits", [])
        for hit in hits:
            if hit.get("approvedSymbol", "").upper() == gene_symbol.upper():
                return hit["id"]
        # Fallback: first target hit
        for hit in hits:
            if hit.get("entity") == "target":
                return hit["id"]
        return None

    def fetch(self, gene: str = "", limit: int = 50, **kwargs) -> Generator[dict[str, Any], None, None]:
        """Fetch target associations and drug data for a gene.

        Args:
            gene: Gene symbol (e.g., HCRT, PSEN1, HCRTR1)
            limit: Max associations/drugs to return
        """
        if not gene:
            return

        ensembl_id = self._resolve_ensembl_id(gene)
        if not ensembl_id:
            return

        data = self._graphql(_TARGET_ASSOCIATIONS_QUERY, {
            "ensemblId": ensembl_id,
            "size": limit,
        })

        target = data.get("target")
        if not target:
            return

        # Yield disease associations
        for assoc in target.get("associatedDiseases", {}).get("rows", []):
            yield {
                "record_type": "disease_association",
                "target": target,
                "association": assoc,
            }

        # Yield known drugs
        for drug_row in target.get("knownDrugs", {}).get("rows", []):
            yield {
                "record_type": "known_drug",
                "target": target,
                "drug_entry": drug_row,
            }

        # Yield tractability summary
        if target.get("tractability"):
            yield {
                "record_type": "tractability",
                "target": target,
                "tractability": target["tractability"],
            }

    def normalize(self, raw_record: dict[str, Any]) -> NormalizedRecord:
        """Normalize OpenTargets records to common schema."""
        record_type = raw_record["record_type"]
        target = raw_record["target"]
        gene = target.get("approvedSymbol", "")

        if record_type == "disease_association":
            assoc = raw_record["association"]
            disease = assoc["disease"]
            score = assoc.get("score", 0)
            return NormalizedRecord(
                record_id=f"{target['id']}_{disease['id']}",
                source_key="opentargets",
                gene_symbol=gene,
                disease=disease.get("name"),
                title=f"{gene} ↔ {disease.get('name')} (score: {score:.3f})",
                summary=f"OpenTargets association score {score:.3f}. "
                        f"Therapeutic areas: {', '.join(ta['name'] for ta in disease.get('therapeuticAreas', []))}",
                payload=assoc,
            )

        elif record_type == "known_drug":
            drug_entry = raw_record["drug_entry"]
            drug = drug_entry.get("drug", {})
            disease = drug_entry.get("disease", {})
            return NormalizedRecord(
                record_id=f"{target['id']}_{drug.get('id', 'unknown')}_{disease.get('id', '')}",
                source_key="opentargets",
                gene_symbol=gene,
                disease=disease.get("name"),
                drug=drug.get("name"),
                title=f"{drug.get('name')} → {gene} (Phase {drug_entry.get('phase', '?')})",
                summary=f"Drug type: {drug.get('drugType')}. "
                        f"Max phase: {drug.get('maximumClinicalTrialPhase')}. "
                        f"Status: {drug_entry.get('status', 'unknown')}.",
                payload=drug_entry,
            )

        elif record_type == "tractability":
            tractability = raw_record["tractability"]
            modalities = {t.get("modality"): t.get("value") for t in tractability}
            return NormalizedRecord(
                record_id=f"{target['id']}_tractability",
                source_key="opentargets",
                gene_symbol=gene,
                title=f"{gene} Druggability/Tractability",
                summary=f"Tractability: {json.dumps(modalities)}",
                payload={"tractability": tractability, "target_id": target["id"]},
            )

        else:
            return NormalizedRecord(
                record_id=f"{target['id']}_{record_type}",
                source_key="opentargets",
                gene_symbol=gene,
                payload=raw_record,
            )

"""cBioPortal Ingestor — Cancer mutation frequencies and molecular profiles.

API: REST at https://www.cbioportal.org/api/
"""
from __future__ import annotations
import json
from typing import Any, Generator
from ..base_ingestor import BaseIngestor, NormalizedRecord
from ..source_registry import SourceConfig


class CbioportalIngestor(BaseIngestor):
    BASE = "https://www.cbioportal.org/api"

    def _apply_auth(self, secret_value: str):
        pass

    def fetch(self, gene: str = "", limit: int = 50, **kwargs) -> Generator[dict[str, Any], None, None]:
        if not gene:
            return
        # Get mutations across all studies
        url = f"{self.BASE}/genes/{gene.upper()}/fetch"
        try:
            resp = self.post(f"{self.BASE}/molecular-profiles/mutations/fetch",
                json={"geneIds": [gene.upper()], "sampleListId": "all"},
                headers={"Content-Type": "application/json"})
        except Exception:
            # Fallback: gene-level info
            resp = self.get(f"{self.BASE}/genes/{gene.upper()}")
            if resp.status_code == 200:
                yield {"record_type": "gene_info", "data": resp.json(), "gene_symbol": gene.upper()}
            return
        if resp.status_code == 200:
            for mut in resp.json()[:limit]:
                yield {"record_type": "mutation", "data": mut, "gene_symbol": gene.upper()}

    def normalize(self, raw: dict[str, Any]) -> NormalizedRecord:
        symbol = raw["gene_symbol"]
        data = raw["data"]
        if raw["record_type"] == "gene_info":
            return NormalizedRecord(
                record_id=f"cbio_gene_{data.get('entrezGeneId', symbol)}",
                source_key="cbioportal",
                gene_symbol=symbol,
                title=f"{symbol} — cBioPortal Gene Entry",
                summary=f"Entrez: {data.get('entrezGeneId')}. Type: {data.get('type', 'unknown')}",
                payload=data,
            )
        else:
            mut_type = data.get("mutationType", "unknown")
            protein = data.get("proteinChange", "")
            study = data.get("studyId", "")
            return NormalizedRecord(
                record_id=f"cbio_{data.get('uniqueMutationId', symbol + '_' + protein)}",
                source_key="cbioportal",
                gene_symbol=symbol,
                disease=data.get("cancerType", None),
                title=f"{symbol} {protein} ({mut_type})",
                summary=f"Study: {study}. Cancer: {data.get('cancerType', '')}. "
                        f"Validation: {data.get('validationStatus', '')}. Center: {data.get('center', '')}",
                payload=data,
            )

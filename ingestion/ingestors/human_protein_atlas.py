"""Human Protein Atlas Ingestor — Tissue expression and protein localization.

API: https://www.proteinatlas.org/{gene}.json
"""
from __future__ import annotations
import json
from typing import Any, Generator
from ..base_ingestor import BaseIngestor, NormalizedRecord
from ..source_registry import SourceConfig


class HumanProteinAtlasIngestor(BaseIngestor):
    BASE = "https://www.proteinatlas.org"

    def _apply_auth(self, secret_value: str):
        pass

    def fetch(self, gene: str = "", limit: int = 50, **kwargs) -> Generator[dict[str, Any], None, None]:
        if not gene:
            return
        resp = self.get(f"{self.BASE}/{gene.upper()}.json")
        data = resp.json()
        if not data:
            return
        # Gene summary
        yield {"record_type": "summary", "data": data, "gene_symbol": gene.upper()}
        # Tissue expression
        for tissue in (data.get("Tissue expression") or [])[:limit]:
            yield {"record_type": "tissue", "data": tissue, "gene_data": data, "gene_symbol": gene.upper()}
        # Subcellular location
        for loc in (data.get("Subcellular location") or [])[:limit]:
            yield {"record_type": "subcellular", "data": loc, "gene_data": data, "gene_symbol": gene.upper()}

    def normalize(self, raw: dict[str, Any]) -> NormalizedRecord:
        symbol = raw["gene_symbol"]
        if raw["record_type"] == "summary":
            d = raw["data"]
            return NormalizedRecord(
                record_id=f"hpa_{symbol}_summary",
                source_key="human_protein_atlas",
                gene_symbol=symbol,
                title=f"{symbol} — Human Protein Atlas",
                summary=f"Name: {d.get('Gene description', '')}. "
                        f"Protein class: {', '.join(d.get('Protein class', []))}. "
                        f"Chromosome: {d.get('Chromosome', '')}. "
                        f"Antibody: {d.get('Antibody', 'N/A')}",
                payload={k: v for k, v in d.items() if k in (
                    "Gene", "Gene description", "Protein class", "Chromosome",
                    "Molecular function", "Biological process", "Disease involvement"
                )},
            )
        elif raw["record_type"] == "tissue":
            t = raw["data"]
            tissue_name = t.get("Tissue", "unknown")
            level = t.get("Level", "")
            return NormalizedRecord(
                record_id=f"hpa_{symbol}_tissue_{tissue_name}",
                source_key="human_protein_atlas",
                gene_symbol=symbol,
                title=f"{symbol} in {tissue_name}: {level}",
                summary=f"Tissue: {tissue_name}. Expression level: {level}. "
                        f"Reliability: {t.get('Reliability', 'N/A')}",
                payload=t,
            )
        else:  # subcellular
            loc = raw["data"]
            location = loc.get("Location", "unknown")
            return NormalizedRecord(
                record_id=f"hpa_{symbol}_loc_{location}",
                source_key="human_protein_atlas",
                gene_symbol=symbol,
                title=f"{symbol} localization: {location}",
                summary=f"Subcellular location: {location}. Reliability: {loc.get('Reliability', 'N/A')}",
                payload=loc,
            )

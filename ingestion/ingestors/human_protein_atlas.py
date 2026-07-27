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

    def _resolve_ensembl_id(self, gene_symbol: str) -> str | None:
        """Resolve gene symbol to Ensembl gene ID via HPA search API."""
        try:
            resp = self.get(
                f"{self.BASE}/api/search_download.php",
                params={"search": gene_symbol, "format": "json", "columns": "g,eg", "compress": "no"},
            )
            results = resp.json()
            for entry in results:
                if entry.get("Gene", "").upper() == gene_symbol.upper():
                    return entry.get("Ensembl")
            return results[0].get("Ensembl") if results else None
        except Exception:
            return None

    def fetch(self, gene: str = "", limit: int = 50, **kwargs) -> Generator[dict[str, Any], None, None]:
        if not gene:
            return
        # Resolve gene symbol to Ensembl ID (HPA URLs require it)
        ensembl_id = self._resolve_ensembl_id(gene.upper())
        if not ensembl_id:
            return
        resp = self.get(f"{self.BASE}/{ensembl_id}.json")
        data = resp.json()
        if not data:
            return
        # Gene summary
        yield {"record_type": "summary", "data": data, "gene_symbol": gene.upper()}
        # Tissue expression (new format: RNA tissue specific nTPM is {tissue: value})
        tissue_ntpm = data.get("RNA tissue specific nTPM") or {}
        if isinstance(tissue_ntpm, dict):
            for tissue_name, ntpm_val in list(tissue_ntpm.items())[:limit]:
                yield {"record_type": "tissue", "data": {"tissue": tissue_name, "nTPM": ntpm_val},
                       "gene_data": data, "gene_symbol": gene.upper()}
        # Subcellular location (new format: flat list of location strings)
        subcell_locs = data.get("Subcellular location") or []
        if isinstance(subcell_locs, list):
            for loc in subcell_locs[:limit]:
                loc_entry = loc if isinstance(loc, dict) else {"location": loc}
                yield {"record_type": "subcellular", "data": loc_entry, "gene_data": data, "gene_symbol": gene.upper()}

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
            tissue_name = t.get("tissue", "unknown")
            ntpm = t.get("nTPM", "")
            return NormalizedRecord(
                record_id=f"hpa_{symbol}_tissue_{tissue_name}",
                source_key="human_protein_atlas",
                gene_symbol=symbol,
                title=f"{symbol} in {tissue_name}: {ntpm} nTPM",
                summary=f"Tissue: {tissue_name}. RNA expression: {ntpm} nTPM. "
                        f"Specificity: {raw.get('gene_data', {}).get('RNA tissue specificity', 'N/A')}",
                payload=t,
            )
        else:  # subcellular
            loc = raw["data"]
            location = loc.get("location", loc.get("Location", "unknown"))
            return NormalizedRecord(
                record_id=f"hpa_{symbol}_loc_{location}",
                source_key="human_protein_atlas",
                gene_symbol=symbol,
                title=f"{symbol} localization: {location}",
                summary=f"Subcellular location: {location}. "
                        f"Main locations: {', '.join(raw.get('gene_data', {}).get('Subcellular main location') or [])}",
                payload=loc,
            )

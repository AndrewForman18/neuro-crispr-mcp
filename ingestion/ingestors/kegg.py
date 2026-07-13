"""KEGG Ingestor — Metabolic and signaling pathways.

API: REST at https://rest.kegg.jp/
"""
from __future__ import annotations
import json
from typing import Any, Generator
from ..base_ingestor import BaseIngestor, NormalizedRecord
from ..source_registry import SourceConfig


class KeggIngestor(BaseIngestor):
    BASE = "https://rest.kegg.jp"

    def _apply_auth(self, secret_value: str):
        pass

    def _parse_kegg_flat(self, text: str) -> dict[str, str]:
        """Parse KEGG flat-file response into key-value pairs."""
        result = {}
        current_key = ""
        for line in text.strip().split("\n"):
            if line and line[0] != " " and line[0] != "\t":
                parts = line.split(None, 1)
                if len(parts) == 2:
                    current_key = parts[0]
                    result[current_key] = parts[1]
                elif len(parts) == 1:
                    current_key = parts[0]
                    result[current_key] = ""
            elif current_key:
                result[current_key] = result[current_key] + " " + line.strip()
        return result

    def fetch(self, gene: str = "", organism: str = "hsa", limit: int = 20, **kwargs) -> Generator[dict[str, Any], None, None]:
        if not gene:
            return
        # Find KEGG gene entry
        resp = self.get(f"{self.BASE}/find/genes/{gene}")
        if resp.status_code != 200 or not resp.text.strip():
            return
        # Parse results (tab-separated: gene_id \t description)
        lines = resp.text.strip().split("\n")
        kegg_gene_id = None
        for line in lines:
            parts = line.split("\t")
            if len(parts) >= 2 and parts[0].startswith(f"{organism}:"):
                kegg_gene_id = parts[0]
                break
        if not kegg_gene_id and lines:
            kegg_gene_id = lines[0].split("\t")[0]
        if not kegg_gene_id:
            return

        # Get pathways for this gene
        resp = self.get(f"{self.BASE}/link/pathway/{kegg_gene_id}")
        if resp.status_code != 200 or not resp.text.strip():
            return
        pathway_ids = []
        for line in resp.text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 2:
                pathway_ids.append(parts[1])

        # Get pathway details
        for pw_id in pathway_ids[:limit]:
            try:
                resp = self.get(f"{self.BASE}/get/{pw_id}")
                if resp.status_code == 200:
                    pw_data = self._parse_kegg_flat(resp.text)
                    yield {"record_type": "pathway", "pathway_id": pw_id,
                           "pathway_data": pw_data, "gene_symbol": gene.upper(),
                           "kegg_gene_id": kegg_gene_id}
            except Exception:
                continue

    def normalize(self, raw: dict[str, Any]) -> NormalizedRecord:
        symbol = raw["gene_symbol"]
        pw_id = raw["pathway_id"]
        pw_data = raw["pathway_data"]
        pw_name = pw_data.get("NAME", "unknown").rstrip(" - Homo sapiens (human)")
        pw_class = pw_data.get("CLASS", "")
        description = pw_data.get("DESCRIPTION", "")

        return NormalizedRecord(
            record_id=f"kegg_{pw_id}_{symbol}",
            source_key="kegg",
            gene_symbol=symbol,
            title=f"{symbol} in {pw_name} ({pw_id})",
            summary=f"Pathway: {pw_name}. Class: {pw_class}. "
                    f"Description: {description[:200]}" if description else
                    f"Pathway: {pw_name}. Class: {pw_class}.",
            payload={"pathway_id": pw_id, "pathway_name": pw_name,
                     "class": pw_class, "kegg_gene_id": raw["kegg_gene_id"]},
        )

"""Reactome Ingestor — Biological pathways and gene membership.

API: REST at https://reactome.org/ContentService/
Note: The /data/pathways/low/entity/ endpoint requires UniProt accessions,
not gene symbols. We resolve gene → UniProt via the UniProt REST API first.
"""
from __future__ import annotations
import json
from typing import Any, Generator
from ..base_ingestor import BaseIngestor, NormalizedRecord
from ..source_registry import SourceConfig


class ReactomeIngestor(BaseIngestor):
    BASE = "https://reactome.org/ContentService"

    def _apply_auth(self, secret_value: str):
        pass

    def _resolve_uniprot(self, gene_symbol: str) -> str | None:
        """Resolve gene symbol to UniProt accession via UniProt REST API."""
        try:
            resp = self.get(
                "https://rest.uniprot.org/uniprotkb/search",
                params={
                    "query": f"gene_exact:{gene_symbol} AND organism_id:9606 AND reviewed:true",
                    "fields": "accession",
                    "format": "json",
                    "size": "1",
                },
            )
            results = resp.json().get("results", [])
            return results[0]["primaryAccession"] if results else None
        except Exception:
            return None

    def fetch(self, gene: str = "", species: str = "Homo sapiens", limit: int = 30, **kwargs) -> Generator[dict[str, Any], None, None]:
        if not gene:
            return
        # Resolve gene symbol to UniProt accession
        accession = self._resolve_uniprot(gene.upper())
        if not accession:
            return
        # Query Reactome mapping endpoint (the entity endpoint no longer accepts symbols)
        resp = self.get(f"{self.BASE}/data/mapping/UniProt/{accession}/pathways")
        if resp.status_code != 200:
            return
        pathways = resp.json()
        for pathway in pathways[:limit]:
            yield {"record_type": "pathway", "pathway": pathway, "gene_symbol": gene.upper()}

    def normalize(self, raw: dict[str, Any]) -> NormalizedRecord:
        symbol = raw["gene_symbol"]
        pw = raw["pathway"]
        pw_id = pw.get("stId", pw.get("dbId", ""))
        pw_name = pw.get("displayName", pw.get("name", "unknown"))
        species = pw.get("speciesName", "")

        return NormalizedRecord(
            record_id=f"reactome_{pw_id}_{symbol}",
            source_key="reactome",
            gene_symbol=symbol,
            title=f"{symbol} in {pw_name}",
            summary=f"Pathway: {pw_name} ({pw_id}). Species: {species}. "
                    f"Compartment: {', '.join(c.get('displayName', '') for c in pw.get('compartment', [])[:3])}. "
                    f"Category: {pw.get('schemaClass', 'Pathway')}",
            payload={"pathway_id": pw_id, "pathway_name": pw_name,
                     "species": species, "schema_class": pw.get("schemaClass", "")},
        )

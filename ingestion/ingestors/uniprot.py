"""UniProt Ingestor — Protein function, structure, and domains.

API: REST at https://rest.uniprot.org/uniprotkb/
"""
from __future__ import annotations
import json
from typing import Any, Generator
from ..base_ingestor import BaseIngestor, NormalizedRecord
from ..source_registry import SourceConfig


class UniprotIngestor(BaseIngestor):
    BASE = "https://rest.uniprot.org/uniprotkb"

    def _apply_auth(self, secret_value: str):
        pass

    def fetch(self, gene: str = "", organism: str = "human", limit: int = 5, **kwargs) -> Generator[dict[str, Any], None, None]:
        if not gene:
            return
        organism_id = "9606" if organism == "human" else "10090"
        query = f"gene_exact:{gene} AND organism_id:{organism_id} AND reviewed:true"
        resp = self.get(f"{self.BASE}/search",
                       params={"query": query, "format": "json", "size": limit,
                               "fields": "accession,id,gene_names,protein_name,organism_name,length,"
                                         "cc_function,cc_subcellular_location,cc_disease,ft_domain,"
                                         "cc_interaction,go_p,go_f,go_c,xref_pdb"})
        data = resp.json()
        for entry in data.get("results", [])[:limit]:
            yield {"record_type": "protein", "data": entry, "gene_symbol": gene.upper()}

    def normalize(self, raw: dict[str, Any]) -> NormalizedRecord:
        symbol = raw["gene_symbol"]
        entry = raw["data"]
        accession = entry.get("primaryAccession", "")
        protein_name = ""
        if entry.get("proteinDescription", {}).get("recommendedName"):
            protein_name = entry["proteinDescription"]["recommendedName"].get("fullName", {}).get("value", "")

        # Extract function
        functions = []
        for comment in entry.get("comments", []):
            if comment.get("commentType") == "FUNCTION":
                for txt in comment.get("texts", []):
                    functions.append(txt.get("value", ""))

        # Extract subcellular location
        locations = []
        for comment in entry.get("comments", []):
            if comment.get("commentType") == "SUBCELLULAR LOCATION":
                for sub in comment.get("subcellularLocations", []):
                    loc = sub.get("location", {}).get("value", "")
                    if loc:
                        locations.append(loc)

        # Extract diseases
        diseases = []
        for comment in entry.get("comments", []):
            if comment.get("commentType") == "DISEASE":
                disease = comment.get("disease", {})
                if disease.get("diseaseId"):
                    diseases.append(disease["diseaseId"])

        return NormalizedRecord(
            record_id=f"uniprot_{accession}",
            source_key="uniprot",
            gene_symbol=symbol,
            disease="; ".join(diseases[:3]) if diseases else None,
            title=f"{symbol} ({accession}) — {protein_name}",
            summary=f"Function: {(functions[0][:300] if functions else 'N/A')}. "
                    f"Location: {', '.join(locations[:5])}. "
                    f"Diseases: {', '.join(diseases[:5]) or 'none'}. "
                    f"Length: {entry.get('sequence', {}).get('length', '?')} aa",
            payload={
                "accession": accession,
                "protein_name": protein_name,
                "function": functions[:2],
                "locations": locations,
                "diseases": diseases,
                "length": entry.get("sequence", {}).get("length"),
                "organism": entry.get("organism", {}).get("scientificName", ""),
            },
        )

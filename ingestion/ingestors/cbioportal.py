"""cBioPortal Ingestor — Somatic mutation profiles from brain/NDD cancer studies.

cBioPortal REST API v3:
  - Resolves HUGO gene symbols to Entrez Gene IDs via /api/genes/{symbol}
  - Uses /api/studies/{studyId}/molecular-profiles to find MUTATION_EXTENDED profiles
  - Fetches mutations via POST .../mutations/fetch with the <studyId>_all sample list
  - Yields per-study mutation summaries + a gene metadata record per gene

API: https://www.cbioportal.org/api/
Docs: https://docs.cbioportal.org/web-api-and-clients/
"""
from __future__ import annotations

import json
import logging
from typing import Any, Generator

from ..base_ingestor import BaseIngestor, NormalizedRecord
from ..source_registry import SourceConfig

logger = logging.getLogger("neuroplex.cbioportal")

# Curated brain/NDD-relevant public studies.  Standard cBioPortal study IDs.
_BRAIN_STUDIES = [
    "gbm_tcga",           # TCGA Glioblastoma Multiforme
    "lgggbm_tcga_pub",    # TCGA LGG + GBM combined
    "brain_cptac_2020",   # CPTAC Brain Cancer 2020
    "mbl_pcgp",           # Pediatric Medulloblastoma (PCGP)
    "gbm_columbia_2019",  # Columbia GBM cohort
    "nbl_amc_2012",       # Neuroblastoma (AMC 2012)
]


class CBioPortalIngestor(BaseIngestor):
    """Ingest somatic mutation data from cBioPortal for brain/NDD cancer studies."""

    def __init__(self, config: SourceConfig):
        super().__init__(config)
        self._gene_cache: dict[str, int] = {}             # symbol → entrezGeneId
        self._profile_cache: dict[str, str | None] = {}   # studyId → mutation profileId

    @property
    def _base(self) -> str:
        return self.config.base_url.rstrip("/")

    def _apply_auth(self, secret_value: str):
        """cBioPortal public API — no authentication required."""
        pass

    # ── Helpers ────────────────────────────────────────────────────────────────

    def _resolve_gene(self, symbol: str) -> dict | None:
        """Return gene metadata dict from cBioPortal (entrezGeneId, chromosome, etc.)."""
        if symbol in self._gene_cache:
            return {"entrezGeneId": self._gene_cache[symbol], "hugoGeneSymbol": symbol}
        try:
            resp = self.get(f"{self._base}/genes/{symbol}",
                            params={"projection": "SUMMARY"})
            data = resp.json()
            if isinstance(data, dict) and "entrezGeneId" in data:
                self._gene_cache[symbol] = data["entrezGeneId"]
                return data
        except Exception as e:
            logger.warning(f"Gene resolution failed for {symbol}: {e}")
        return None

    def _get_mutation_profile(self, study_id: str) -> str | None:
        """Return the MUTATION_EXTENDED molecular profile ID for a study (cached).

        Uses /api/studies/{studyId}/molecular-profiles — the study-scoped endpoint.
        The global /api/molecular-profiles endpoint is not filtered by studyId.
        """
        if study_id in self._profile_cache:
            return self._profile_cache[study_id]
        try:
            resp = self.get(f"{self._base}/studies/{study_id}/molecular-profiles",
                            params={"projection": "SUMMARY"})
            for p in resp.json():
                if p.get("molecularAlterationType") == "MUTATION_EXTENDED":
                    self._profile_cache[study_id] = p["molecularProfileId"]
                    return p["molecularProfileId"]
        except Exception as e:
            logger.debug(f"Profile lookup skipped for {study_id}: {e}")
        self._profile_cache[study_id] = None
        return None

    def _fetch_mutations(self, profile_id: str, study_id: str,
                         entrez_id: int) -> list[dict]:
        """Fetch all mutations for one gene in one molecular profile.

        POST .../mutations/fetch with the <studyId>_all sample list retrieves
        every mutation in the cohort for this gene without enumerating samples.
        """
        try:
            resp = self.post(
                f"{self._base}/molecular-profiles/{profile_id}/mutations/fetch",
                params={"projection": "SUMMARY", "pageSize": 500},
                json={
                    "sampleListId": f"{study_id}_all",
                    "entrezGeneIds": [entrez_id],
                },
            )
            result = resp.json()
            return result if isinstance(result, list) else []
        except Exception as e:
            logger.debug(f"Mutation fetch skipped — {profile_id}/{entrez_id}: {e}")
            return []

    # ── BaseIngestor interface ──────────────────────────────────────────────────

    def fetch(self, gene: str = "", **kwargs) -> Generator[dict[str, Any], None, None]:
        """Fetch mutation data for a gene across brain/NDD cancer studies.

        Args:
            gene: HUGO gene symbol (e.g. HCRT, PSEN1, APP)
        Yields:
            dicts with keys: record_type, gene_symbol, gene_meta,
                             study_id, profile_id, mutations
        """
        if not gene:
            return

        gene_meta = self._resolve_gene(gene)
        if not gene_meta:
            logger.warning(f"cBioPortal: no gene entry for {gene}")
            return

        entrez_id: int = gene_meta["entrezGeneId"]

        # Always emit a gene_info record so the table has an entry per gene
        yield {
            "record_type": "gene_info",
            "gene_symbol": gene,
            "gene_meta": gene_meta,
            "study_id": None,
            "profile_id": None,
            "mutations": [],
        }

        for study_id in _BRAIN_STUDIES:
            profile_id = self._get_mutation_profile(study_id)
            if not profile_id:
                continue
            mutations = self._fetch_mutations(profile_id, study_id, entrez_id)
            if mutations:
                yield {
                    "record_type": "study_mutations",
                    "gene_symbol": gene,
                    "gene_meta": gene_meta,
                    "study_id": study_id,
                    "profile_id": profile_id,
                    "mutations": mutations,
                }

    def normalize(self, raw: dict[str, Any]) -> NormalizedRecord:
        """Normalise a cBioPortal fetch record into the common NeuroPlex schema."""
        gene      = raw["gene_symbol"]
        gene_meta = raw.get("gene_meta") or {}
        rtype     = raw["record_type"]

        if rtype == "gene_info":
            return NormalizedRecord(
                record_id  = f"cbio_gene_{gene_meta.get('entrezGeneId', gene)}",
                source_key = "cbioportal",
                gene_symbol= gene,
                title      = f"{gene} — cBioPortal Gene Entry",
                summary    = (
                    f"Entrez: {gene_meta.get('entrezGeneId')} | "
                    f"Chr: {gene_meta.get('chromosome', '?')} "
                    f"{gene_meta.get('cytoband', '')} | "
                    f"Type: {gene_meta.get('type', 'unknown')}"
                ),
                payload = gene_meta,
            )

        # record_type == "study_mutations"
        study_id  = raw["study_id"]
        mutations = raw["mutations"]

        mut_types: dict[str, int] = {}
        samples: set[str] = set()
        protein_changes: list[str] = []
        for m in mutations:
            mt = m.get("mutationType") or m.get("mutationStatus") or "Unknown"
            mut_types[mt] = mut_types.get(mt, 0) + 1
            if sid := m.get("sampleId"):
                samples.add(sid)
            if pc := m.get("proteinChange"):
                protein_changes.append(pc)

        n_mut    = len(mutations)
        n_samp   = len(samples)
        dominant = max(mut_types, key=mut_types.get) if mut_types else "Unknown"

        return NormalizedRecord(
            record_id  = f"cbio_{gene}_{study_id}",
            source_key = "cbioportal",
            gene_symbol= gene,
            title      = f"{gene} somatic mutations in {study_id} — {n_mut} events, {n_samp} samples",
            summary    = (
                f"Gene: {gene} | Study: {study_id} | "
                f"Mutations: {n_mut} | Samples: {n_samp} | "
                f"Dominant type: {dominant} | "
                f"Types: {json.dumps(mut_types)} | "
                f"Top protein changes: {', '.join(protein_changes[:10])}"
            ),
            payload = {
                "gene_symbol":         gene,
                "entrez_id":           gene_meta.get("entrezGeneId"),
                "study_id":            study_id,
                "molecular_profile_id": raw["profile_id"],
                "n_mutations":         n_mut,
                "n_samples_affected":  n_samp,
                "mutation_types":      mut_types,
                "top_protein_changes": protein_changes[:50],
                "mutations_sample":    mutations[:30],
            },
        )

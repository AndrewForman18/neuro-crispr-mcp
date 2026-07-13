"""gnomAD Ingestor — Population allele frequencies and gene constraint.

gnomAD (Genome Aggregation Database) v4.1:
- Population allele frequencies across 800K+ exomes/genomes
- Gene constraint metrics (pLI, LOEUF, mis_z) for loss-of-function intolerance
- Variant-level annotations (consequence, clinical significance flags)
- Population-specific frequencies (NFE, EAS, AFR, AMR, SAS, ASJ, FIN, MID)

API: GraphQL at https://gnomad.broadinstitute.org/api
Docs: https://gnomad.broadinstitute.org/help

Key for NeuroPlex:
- Gene constraint scores identify intolerant-to-LoF targets (NDD/orexin pathway)
- Population frequencies contextualize ClinVar pathogenic variants
- Variant impact assessment for CRISPR knockout phenotype interpretation
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from typing import Any, Generator

from ..base_ingestor import BaseIngestor, NormalizedRecord
from ..source_registry import SourceConfig


# ═══════════════════════════════════════════════════════════════════════════════
# GraphQL Queries
# ═══════════════════════════════════════════════════════════════════════════════

_GENE_CONSTRAINT_QUERY = """
query GeneConstraint($geneSymbol: String!, $referenceGenome: ReferenceGenomeId!) {
  gene(gene_symbol: $geneSymbol, reference_genome: $referenceGenome) {
    gene_id
    symbol
    name
    chrom
    start
    stop
    strand
    gnomad_constraint {
      exp_lof
      exp_mis
      exp_syn
      obs_lof
      obs_mis
      obs_syn
      oe_lof
      oe_lof_lower
      oe_lof_upper
      oe_mis
      oe_mis_lower
      oe_mis_upper
      oe_syn
      oe_syn_lower
      oe_syn_upper
      lof_z
      mis_z
      syn_z
      pLI
      flags
    }
    clinvar_variants {
      clinical_significance
      clinvar_variation_id
      gnomad {
        exome {
          ac
          an
        }
        genome {
          ac
          an
        }
      }
      gold_stars
      hgvsc
      hgvsp
      in_gnomad
      major_consequence
      pos
      review_status
      variant_id
    }
  }
}
"""

_GENE_VARIANTS_QUERY = """
query GeneVariants($geneId: String!, $datasetId: DatasetId!, $referenceGenome: ReferenceGenomeId!) {
  gene(gene_id: $geneId, reference_genome: $referenceGenome) {
    variants(dataset: $datasetId) {
      variant_id
      pos
      consequence
      hgvsc
      hgvsp
      lof
      lof_filter
      lof_flags
      exome {
        ac
        an
        af
        homozygote_count
        populations {
          id
          ac
          an
          af
        }
      }
      genome {
        ac
        an
        af
        homozygote_count
        populations {
          id
          ac
          an
          af
        }
      }
    }
  }
}
"""


# ═══════════════════════════════════════════════════════════════════════════════
# Ingestor
# ═══════════════════════════════════════════════════════════════════════════════

class GnomadIngestor(BaseIngestor):
    """Ingest gene constraint and variant frequency data from gnomAD.

    Records produced per gene:
    1. Gene constraint record (pLI, LOEUF, o/e ratios)
    2. ClinVar variant records (pathogenic variants with gnomAD frequencies)
    3. Loss-of-function variant records (pLoF variants with population AFs)
    """

    GNOMAD_API = "https://gnomad.broadinstitute.org/api"
    REFERENCE_GENOME = "GRCh38"
    DATASET_ID = "gnomad_r4"  # gnomAD v4.1

    def _apply_auth(self, secret_value: str):
        """gnomAD API is open — no auth required."""
        pass

    def _graphql(self, query: str, variables: dict) -> dict:
        """Execute a GraphQL query against gnomAD API."""
        resp = self.post(
            self.GNOMAD_API,
            json={"query": query, "variables": variables},
            headers={"Content-Type": "application/json"},
        )
        data = resp.json()
        if "errors" in data:
            errors = data["errors"]
            msg = "; ".join(e.get("message", str(e)) for e in errors)
            raise RuntimeError(f"gnomAD GraphQL error: {msg}")
        return data.get("data", {})

    def fetch(self, gene: str = "", include_variants: bool = False, limit: int = 50, **kwargs) -> Generator[dict[str, Any], None, None]:
        """Fetch gene constraint + ClinVar + optionally pLoF variants.

        Args:
            gene: Gene symbol (e.g., HCRT, PSEN1, HCRTR1)
            include_variants: If True, also fetch individual pLoF variants (slower)
            limit: Max variants to return when include_variants=True
        """
        if not gene:
            return

        # 1. Gene constraint + ClinVar variants
        data = self._graphql(_GENE_CONSTRAINT_QUERY, {
            "geneSymbol": gene.upper(),
            "referenceGenome": self.REFERENCE_GENOME,
        })

        gene_data = data.get("gene")
        if not gene_data:
            return

        # Yield constraint record
        constraint = gene_data.get("gnomad_constraint")
        if constraint:
            yield {
                "record_type": "gene_constraint",
                "gene": gene_data,
                "constraint": constraint,
            }

        # Yield ClinVar variant records
        clinvar_variants = gene_data.get("clinvar_variants") or []
        for variant in clinvar_variants[:limit]:
            yield {
                "record_type": "clinvar_variant",
                "gene": gene_data,
                "variant": variant,
            }

        # 2. Optionally fetch pLoF variants
        if include_variants and gene_data.get("gene_id"):
            try:
                var_data = self._graphql(_GENE_VARIANTS_QUERY, {
                    "geneId": gene_data["gene_id"],
                    "datasetId": self.DATASET_ID,
                    "referenceGenome": self.REFERENCE_GENOME,
                })
                variants = var_data.get("gene", {}).get("variants") or []
                # Filter to pLoF (high-confidence loss-of-function)
                plof_variants = [
                    v for v in variants
                    if v.get("lof") == "HC"  # High-confidence LoF
                    or v.get("consequence") in (
                        "frameshift_variant", "stop_gained",
                        "splice_acceptor_variant", "splice_donor_variant",
                    )
                ]
                for variant in plof_variants[:limit]:
                    yield {
                        "record_type": "plof_variant",
                        "gene": gene_data,
                        "variant": variant,
                    }
            except Exception as e:
                # Variant fetch is optional — don't fail the whole gene
                import logging
                logging.getLogger("neuroplex.ingestor").warning(
                    f"Failed to fetch variants for {gene}: {e}"
                )

    def normalize(self, raw_record: dict[str, Any]) -> NormalizedRecord:
        """Normalize gnomAD records to common schema."""
        record_type = raw_record["record_type"]
        gene_data = raw_record["gene"]
        symbol = gene_data.get("symbol", "").upper()
        gene_id = gene_data.get("gene_id", "")

        if record_type == "gene_constraint":
            constraint = raw_record["constraint"]
            pli = constraint.get("pLI", 0)
            loeuf = constraint.get("oe_lof_upper", None)
            mis_z = constraint.get("mis_z", None)

            # Interpret constraint
            if pli and pli > 0.9:
                interpretation = "highly intolerant to LoF (haploinsufficient)"
            elif pli and pli > 0.5:
                interpretation = "moderately constrained"
            else:
                interpretation = "tolerant to LoF"

            return NormalizedRecord(
                record_id=f"{gene_id}_constraint",
                source_key="gnomad",
                gene_symbol=symbol,
                title=f"{symbol} Gene Constraint (gnomAD v4)",
                summary=(
                    f"pLI={pli:.4f} ({interpretation}). "
                    f"LOEUF={loeuf:.3f}. " if loeuf else f"pLI={pli:.4f} ({interpretation}). "
                    f"mis_z={mis_z:.2f}. " if mis_z else ""
                    f"obs/exp LoF: {constraint.get('obs_lof', '?')}/{constraint.get('exp_lof', '?'):.1f}. "
                    f"Location: chr{gene_data.get('chrom', '?')}:{gene_data.get('start', '?')}-{gene_data.get('stop', '?')}"
                ),
                payload={
                    "constraint": constraint,
                    "gene_id": gene_id,
                    "chrom": gene_data.get("chrom"),
                    "start": gene_data.get("start"),
                    "stop": gene_data.get("stop"),
                    "strand": gene_data.get("strand"),
                    "gene_name": gene_data.get("name"),
                },
            )

        elif record_type == "clinvar_variant":
            variant = raw_record["variant"]
            variant_id = variant.get("variant_id", "unknown")
            clin_sig = variant.get("clinical_significance", "unknown")
            consequence = variant.get("major_consequence", "")
            hgvsp = variant.get("hgvsp", "")
            hgvsc = variant.get("hgvsc", "")

            # Extract gnomAD frequency (compute AF from ac/an)
            gnomad = variant.get("gnomad", {})
            exome = gnomad.get("exome") or {}
            genome = gnomad.get("genome") or {}
            exome_af = (exome["ac"] / exome["an"]) if exome.get("ac") and exome.get("an") and exome["an"] > 0 else None
            genome_af = (genome["ac"] / genome["an"]) if genome.get("ac") and genome.get("an") and genome["an"] > 0 else None
            best_af = exome_af or genome_af or 0

            return NormalizedRecord(
                record_id=f"{gene_id}_clinvar_{variant.get('clinvar_variation_id', variant_id)}",
                source_key="gnomad",
                gene_symbol=symbol,
                title=f"{symbol} {hgvsp or hgvsc or variant_id} — {clin_sig}",
                summary=(
                    f"ClinVar: {clin_sig} ({variant.get('review_status', '')}, "
                    f"{variant.get('gold_stars', 0)} stars). "
                    f"Consequence: {consequence}. "
                    f"gnomAD AF: {best_af:.6f}" + (
                        " (rare)" if best_af and best_af < 0.001 else
                        " (ultra-rare)" if best_af and best_af < 0.0001 else
                        " (common)" if best_af and best_af > 0.01 else ""
                    ) + f". In gnomAD: {'Yes' if variant.get('in_gnomad') else 'No'}"
                ),
                payload=variant,
            )

        elif record_type == "plof_variant":
            variant = raw_record["variant"]
            variant_id = variant.get("variant_id", "unknown")
            consequence = variant.get("consequence", "")
            lof_flag = variant.get("lof", "")
            lof_filter = variant.get("lof_filter", "")

            # Get best allele frequency
            exome = variant.get("exome") or {}
            genome = variant.get("genome") or {}
            best_af = exome.get("af") or genome.get("af") or 0
            best_ac = exome.get("ac") or genome.get("ac") or 0
            hom = exome.get("homozygote_count", 0) + genome.get("homozygote_count", 0)

            # Population-specific frequencies
            populations = exome.get("populations") or genome.get("populations") or []
            pop_summary = {p["id"]: p["af"] for p in populations if p.get("af") and p["af"] > 0}

            return NormalizedRecord(
                record_id=f"{gene_id}_plof_{variant_id}",
                source_key="gnomad",
                gene_symbol=symbol,
                title=f"{symbol} pLoF: {variant_id} ({consequence})",
                summary=(
                    f"Predicted LoF ({lof_flag}): {consequence}. "
                    f"AF={best_af:.6f}, AC={best_ac}, Hom={hom}. "
                    f"Filter: {lof_filter or 'PASS'}. "
                    f"Pops: {json.dumps(pop_summary) if pop_summary else 'rare across all'}"
                ),
                payload={
                    "variant_id": variant_id,
                    "consequence": consequence,
                    "lof": lof_flag,
                    "lof_filter": lof_filter,
                    "lof_flags": variant.get("lof_flags"),
                    "hgvsc": variant.get("hgvsc"),
                    "hgvsp": variant.get("hgvsp"),
                    "exome_af": exome.get("af"),
                    "exome_ac": exome.get("ac"),
                    "genome_af": genome.get("af"),
                    "genome_ac": genome.get("ac"),
                    "homozygote_count": hom,
                    "populations": pop_summary,
                },
            )

        else:
            return NormalizedRecord(
                record_id=f"{gene_id}_{record_type}",
                source_key="gnomad",
                gene_symbol=symbol,
                payload=raw_record,
            )

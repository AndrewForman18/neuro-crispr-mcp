"""NeuroPlex Query Module — Unified gene profiling with Plotly graphics.

Provides cross-source gene intelligence queries and chart generators for the
Neuro CRISPR MCP Streamlit UI. Queries all populated neuroplex_* Delta tables
and returns structured data + Plotly figure JSON for rendering.

Sources:
  - OpenTargets: disease associations with evidence scores
  - gnomAD v4: ClinVar variants + gene constraint metrics
  - KEGG: pathway memberships
  - cBioPortal: somatic mutations
  - UniProt: protein function annotations
  - NCBI GTR: genetic testing registry entries
"""

import json
import logging
from typing import Any

import plotly.graph_objects as go
import numpy as np

from config.neuroplex_config import load_config

logger = logging.getLogger(__name__)

# ─── Table Definitions ───────────────────────────────────────────────────────
CFG = load_config()
CATALOG = CFG.catalog
SCHEMA = CFG.schema

TABLE_OPENTARGETS = f"{CATALOG}.{SCHEMA}.neuroplex_opentargets"
TABLE_GNOMAD = f"{CATALOG}.{SCHEMA}.neuroplex_gnomad"
TABLE_KEGG = f"{CATALOG}.{SCHEMA}.neuroplex_kegg"
TABLE_CBIOPORTAL = f"{CATALOG}.{SCHEMA}.neuroplex_cbioportal"
TABLE_UNIPROT = f"{CATALOG}.{SCHEMA}.neuroplex_uniprot"
TABLE_NCBI_GTR = f"{CATALOG}.{SCHEMA}.neuroplex_ncbi_gtr"

# Color palette matching the app's dark theme
COLORS = {
    "primary": "#FF4B4B",
    "secondary": "#1f77b4",
    "accent": "#00cc96",
    "warning": "#ffa15a",
    "muted": "#636efa",
    "bg": "#0e1117",
    "card_bg": "#1e1e2e",
    "text": "#fafafa",
    "grid": "#2a2a3e",
}

DISEASE_PALETTE = [
    "#636efa", "#ef553b", "#00cc96", "#ab63fa", "#ffa15a",
    "#19d3f3", "#ff6692", "#b6e880", "#ff97ff", "#fecb52",
]


# ─── Shared Plotly Layout ────────────────────────────────────────────────────
def _base_layout(**overrides) -> dict:
    """Dark-themed Plotly layout matching the Streamlit app."""
    layout = dict(
        template="plotly_dark",
        paper_bgcolor=COLORS["bg"],
        plot_bgcolor=COLORS["card_bg"],
        font=dict(family="Inter, sans-serif", color=COLORS["text"], size=12),
        margin=dict(l=60, r=30, t=50, b=50),
        xaxis=dict(gridcolor=COLORS["grid"], zeroline=False),
        yaxis=dict(gridcolor=COLORS["grid"], zeroline=False),
    )
    layout.update(overrides)
    return layout


# ─── Query Functions ─────────────────────────────────────────────────────────

def query_neuroplex_gene_profile(gene_symbol: str, _execute_query_fn=None) -> str:
    """Unified NeuroPlex gene profile across all sources.

    Returns aggregated intelligence from OpenTargets, gnomAD, KEGG, cBioPortal,
    UniProt, and NCBI GTR for a single gene. Includes summary statistics and
    Plotly chart JSON for each data facet.

    Args:
        gene_symbol: Gene symbol (UPPER CASE, e.g. HCRT, PSEN1, APP, SOD1).
        _execute_query_fn: Injected SQL execution function (from server.py).
    """
    if not _execute_query_fn:
        return json.dumps({"error": "No SQL executor provided."})

    gene = gene_symbol.upper().strip()
    profile = {"gene": gene, "sources": {}}

    # ── OpenTargets: Disease associations ──
    ot_query = f"""
        SELECT disease,
               payload:score::DOUBLE AS association_score,
               payload:disease.therapeuticAreas[0].name::STRING AS therapeutic_area
        FROM {TABLE_OPENTARGETS}
        WHERE UPPER(gene_symbol) = '{gene}'
        ORDER BY payload:score::DOUBLE DESC NULLS LAST
        LIMIT 20
    """
    try:
        ot_results = _execute_query_fn(ot_query, timeout=30)
        profile["sources"]["opentargets"] = {
            "n_associations": len(ot_results),
            "top_diseases": ot_results[:10],
        }
    except Exception as e:
        logger.warning(f"OpenTargets query failed: {e}")
        profile["sources"]["opentargets"] = {"error": str(e)}

    # ── gnomAD: Variants + constraint ──
    gn_variants_query = f"""
        SELECT title,
               payload:clinical_significance::STRING AS clinical_significance,
               payload:hgvsp::STRING AS protein_change,
               payload:major_consequence::STRING AS consequence,
               payload:pos::INT AS position,
               payload:gnomad.exome.ac::INT AS exome_ac,
               payload:gnomad.exome.an::INT AS exome_an,
               payload:gold_stars::INT AS review_stars
        FROM {TABLE_GNOMAD}
        WHERE UPPER(gene_symbol) = '{gene}'
          AND title NOT LIKE '%Constraint%'
        ORDER BY payload:gold_stars::INT DESC NULLS LAST
        LIMIT 50
    """
    gn_constraint_query = f"""
        SELECT payload:constraint.pLI::DOUBLE AS pLI,
               payload:constraint.oe_lof::DOUBLE AS oe_lof,
               payload:constraint.oe_lof_upper::DOUBLE AS oe_lof_upper,
               payload:constraint.oe_mis::DOUBLE AS oe_mis,
               payload:constraint.mis_z::DOUBLE AS mis_z,
               payload:constraint.lof_z::DOUBLE AS lof_z,
               payload:constraint.obs_lof::INT AS obs_lof,
               payload:constraint.exp_lof::DOUBLE AS exp_lof,
               payload:constraint.obs_mis::INT AS obs_mis,
               payload:constraint.exp_mis::DOUBLE AS exp_mis,
               payload:gene_name::STRING AS full_gene_name
        FROM {TABLE_GNOMAD}
        WHERE UPPER(gene_symbol) = '{gene}'
          AND title LIKE '%Constraint%'
        LIMIT 1
    """
    try:
        gn_variants = _execute_query_fn(gn_variants_query, timeout=30)
        gn_constraint = _execute_query_fn(gn_constraint_query, timeout=30)
        profile["sources"]["gnomad"] = {
            "n_variants": len(gn_variants),
            "constraint": gn_constraint[0] if gn_constraint else None,
            "variants": gn_variants[:20],
        }
    except Exception as e:
        logger.warning(f"gnomAD query failed: {e}")
        profile["sources"]["gnomad"] = {"error": str(e)}

    # ── KEGG: Pathway memberships ──
    kegg_query = f"""
        SELECT title, summary,
               payload:pathway_id::STRING AS pathway_id,
               payload:pathway_name::STRING AS pathway_name
        FROM {TABLE_KEGG}
        WHERE UPPER(gene_symbol) = '{gene}'
        ORDER BY title
    """
    try:
        kegg_results = _execute_query_fn(kegg_query, timeout=30)
        profile["sources"]["kegg"] = {
            "n_pathways": len(kegg_results),
            "pathways": kegg_results,
        }
    except Exception as e:
        logger.warning(f"KEGG query failed: {e}")
        profile["sources"]["kegg"] = {"error": str(e)}

    # ── cBioPortal: Somatic mutations ──
    cbio_query = f"""
        SELECT title, disease, summary,
               payload:mutation_type::STRING AS mutation_type,
               payload:protein_change::STRING AS protein_change,
               payload:cancer_study::STRING AS cancer_study
        FROM {TABLE_CBIOPORTAL}
        WHERE UPPER(gene_symbol) = '{gene}'
        LIMIT 20
    """
    try:
        cbio_results = _execute_query_fn(cbio_query, timeout=30)
        profile["sources"]["cbioportal"] = {
            "n_mutations": len(cbio_results),
            "mutations": cbio_results,
        }
    except Exception as e:
        logger.warning(f"cBioPortal query failed: {e}")
        profile["sources"]["cbioportal"] = {"error": str(e)}

    # ── UniProt: Protein annotations ──
    uniprot_query = f"""
        SELECT title, summary,
               payload:protein_name::STRING AS protein_name,
               payload:function::STRING AS function_desc,
               payload:subcellular_location::STRING AS location
        FROM {TABLE_UNIPROT}
        WHERE UPPER(gene_symbol) = '{gene}'
        LIMIT 5
    """
    try:
        uniprot_results = _execute_query_fn(uniprot_query, timeout=30)
        profile["sources"]["uniprot"] = {
            "n_entries": len(uniprot_results),
            "annotations": uniprot_results,
        }
    except Exception as e:
        logger.warning(f"UniProt query failed: {e}")
        profile["sources"]["uniprot"] = {"error": str(e)}

    # ── NCBI GTR: Genetic tests ──
    gtr_query = f"""
        SELECT title, summary, disease,
               payload:test_type::STRING AS test_type,
               payload:lab_name::STRING AS lab_name
        FROM {TABLE_NCBI_GTR}
        WHERE UPPER(gene_symbol) = '{gene}'
        LIMIT 20
    """
    try:
        gtr_results = _execute_query_fn(gtr_query, timeout=30)
        profile["sources"]["ncbi_gtr"] = {
            "n_tests": len(gtr_results),
            "tests": gtr_results[:10],
        }
    except Exception as e:
        logger.warning(f"NCBI GTR query failed: {e}")
        profile["sources"]["ncbi_gtr"] = {"error": str(e)}

    # ── Summary stats ──
    profile["summary"] = {
        "total_sources_with_data": sum(
            1 for v in profile["sources"].values()
            if "error" not in v and any(
                v.get(k, 0) for k in ["n_associations", "n_variants", "n_pathways",
                                       "n_mutations", "n_entries", "n_tests"]
            )
        ),
        "total_records": sum(
            v.get(k, 0)
            for v in profile["sources"].values()
            if "error" not in v
            for k in ["n_associations", "n_variants", "n_pathways",
                      "n_mutations", "n_entries", "n_tests"]
            if k in v
        ),
    }

    return json.dumps(profile, default=str)


def query_neuroplex_disease_landscape(disease_filter: str = "", limit: int = 15, _execute_query_fn=None) -> str:
    """Query NeuroPlex disease-gene landscape from OpenTargets.

    Returns all gene-disease associations, optionally filtered by disease name.
    Useful for seeing which genes in the panel are linked to a specific condition.

    Args:
        disease_filter: Disease name substring (e.g. narcolepsy, Alzheimer, Parkinson).
        limit: Max results per gene.
        _execute_query_fn: Injected SQL execution function.
    """
    if not _execute_query_fn:
        return json.dumps({"error": "No SQL executor provided."})

    where = "WHERE 1=1"
    if disease_filter:
        escaped = disease_filter.replace("'", "''")
        where += f" AND LOWER(disease) LIKE LOWER('%{escaped}%')"

    query = f"""
        SELECT gene_symbol, disease,
               payload:score::DOUBLE AS score,
               payload:disease.therapeuticAreas[0].name::STRING AS therapeutic_area
        FROM {TABLE_OPENTARGETS}
        {where}
        ORDER BY payload:score::DOUBLE DESC NULLS LAST
        LIMIT {min(limit * 15, 200)}
    """
    try:
        results = _execute_query_fn(query, timeout=30)
        return json.dumps({
            "dataset": "NeuroPlex (OpenTargets)",
            "disease_filter": disease_filter or "all",
            "associations_found": len(results),
            "results": results,
        }, default=str)
    except Exception as e:
        return json.dumps({"error": str(e)})


def query_neuroplex_variant_summary(gene_symbol: str = "", _execute_query_fn=None) -> str:
    """Query NeuroPlex variant summary across all genes from gnomAD.

    Returns clinical significance distribution and constraint metrics.

    Args:
        gene_symbol: Optional gene filter (UPPER CASE). If empty, returns panel-wide summary.
        _execute_query_fn: Injected SQL execution function.
    """
    if not _execute_query_fn:
        return json.dumps({"error": "No SQL executor provided."})

    gene_filter = ""
    if gene_symbol:
        gene_filter = f"AND UPPER(gene_symbol) = '{gene_symbol.upper().strip()}'"

    # Clinical significance distribution
    sig_query = f"""
        SELECT gene_symbol,
               payload:clinical_significance::STRING AS clinical_significance,
               COUNT(*) AS variant_count
        FROM {TABLE_GNOMAD}
        WHERE title NOT LIKE '%Constraint%'
          {gene_filter}
        GROUP BY gene_symbol, payload:clinical_significance::STRING
        ORDER BY gene_symbol, variant_count DESC
    """
    # Constraint metrics for all genes
    constraint_query = f"""
        SELECT gene_symbol,
               payload:constraint.pLI::DOUBLE AS pLI,
               payload:constraint.oe_lof::DOUBLE AS oe_lof,
               payload:constraint.oe_lof_upper::DOUBLE AS oe_lof_upper,
               payload:constraint.mis_z::DOUBLE AS mis_z,
               payload:constraint.lof_z::DOUBLE AS lof_z
        FROM {TABLE_GNOMAD}
        WHERE title LIKE '%Constraint%'
          {gene_filter}
        ORDER BY payload:constraint.pLI::DOUBLE DESC NULLS LAST
    """
    try:
        sig_results = _execute_query_fn(sig_query, timeout=30)
        constraint_results = _execute_query_fn(constraint_query, timeout=30)
        return json.dumps({
            "dataset": "NeuroPlex (gnomAD v4)",
            "gene_filter": gene_symbol or "all panel genes",
            "variant_classification": sig_results,
            "gene_constraint": constraint_results,
        }, default=str)
    except Exception as e:
        return json.dumps({"error": str(e)})


# ─── Chart Generators ────────────────────────────────────────────────────────

def chart_disease_associations(profile_json: str) -> go.Figure | None:
    """Horizontal bar chart of disease association scores from OpenTargets.

    Args:
        profile_json: JSON string from query_neuroplex_gene_profile.
    """
    try:
        profile = json.loads(profile_json)
        ot_data = profile.get("sources", {}).get("opentargets", {})
        diseases = ot_data.get("top_diseases", [])
        if not diseases:
            return None

        gene = profile.get("gene", "")
        names = [d.get("disease", "Unknown")[:40] for d in diseases]
        scores = [float(d.get("association_score") or 0) for d in diseases]
        areas = [d.get("therapeutic_area", "Unknown") for d in diseases]

        # Reverse for horizontal bar (top items at top)
        names.reverse()
        scores.reverse()
        areas.reverse()

        # Color by therapeutic area
        unique_areas = list(dict.fromkeys(areas))
        color_map = {a: DISEASE_PALETTE[i % len(DISEASE_PALETTE)] for i, a in enumerate(unique_areas)}
        colors = [color_map[a] for a in areas]

        fig = go.Figure()
        fig.add_trace(go.Bar(
            x=scores,
            y=names,
            orientation="h",
            marker=dict(color=colors, line=dict(width=0)),
            text=[f"{s:.3f}" for s in scores],
            textposition="outside",
            hovertemplate="<b>%{y}</b><br>Score: %{x:.4f}<extra></extra>",
        ))

        fig.update_layout(
            **_base_layout(
                title=dict(text=f"{gene} — Disease Associations (OpenTargets)", font=dict(size=16)),
                xaxis_title="Association Score",
                yaxis_title=None,
                height=max(350, len(names) * 32 + 100),
                showlegend=False,
            )
        )
        return fig
    except Exception as e:
        logger.warning(f"chart_disease_associations error: {e}")
        return None


def chart_variant_landscape(profile_json: str) -> go.Figure | None:
    """Stacked bar chart of ClinVar variant classification from gnomAD.

    Args:
        profile_json: JSON string from query_neuroplex_gene_profile.
    """
    try:
        profile = json.loads(profile_json)
        gn_data = profile.get("sources", {}).get("gnomad", {})
        variants = gn_data.get("variants", [])
        if not variants:
            return None

        gene = profile.get("gene", "")

        # Count by clinical significance
        sig_counts: dict[str, int] = {}
        for v in variants:
            sig = v.get("clinical_significance") or "Unknown"
            sig_counts[sig] = sig_counts.get(sig, 0) + 1

        # Sort by count descending
        sorted_sigs = sorted(sig_counts.items(), key=lambda x: x[1], reverse=True)
        labels = [s[0] for s in sorted_sigs]
        counts = [s[1] for s in sorted_sigs]

        # Color map for clinical significance
        sig_colors = {
            "Pathogenic": "#ef553b",
            "Likely pathogenic": "#ff6692",
            "Uncertain significance": "#ffa15a",
            "Likely benign": "#00cc96",
            "Benign": "#19d3f3",
            "Unknown": "#636efa",
        }
        colors = [sig_colors.get(l, "#636efa") for l in labels]

        fig = go.Figure()
        fig.add_trace(go.Bar(
            x=labels,
            y=counts,
            marker=dict(color=colors, line=dict(width=1, color=COLORS["grid"])),
            text=counts,
            textposition="outside",
            hovertemplate="<b>%{x}</b><br>Count: %{y}<extra></extra>",
        ))

        fig.update_layout(
            **_base_layout(
                title=dict(text=f"{gene} — ClinVar Variant Classification (gnomAD v4)", font=dict(size=16)),
                xaxis_title="Clinical Significance",
                yaxis_title="Variant Count",
                height=400,
                showlegend=False,
            )
        )
        return fig
    except Exception as e:
        logger.warning(f"chart_variant_landscape error: {e}")
        return None


def chart_gene_constraint(profile_json: str) -> go.Figure | None:
    """Radar chart of gene constraint metrics from gnomAD.

    Visualizes pLI, LOEUF, missense Z, LoF Z, and o/e ratios.

    Args:
        profile_json: JSON string from query_neuroplex_gene_profile.
    """
    try:
        profile = json.loads(profile_json)
        gn_data = profile.get("sources", {}).get("gnomad", {})
        constraint = gn_data.get("constraint")
        if not constraint:
            return None

        gene = profile.get("gene", "")

        # Metrics for radar (normalize to 0-1 range)
        pli = float(constraint.get("pLI") or 0)
        oe_lof = float(constraint.get("oe_lof") or 1)
        oe_mis = float(constraint.get("oe_mis") or 1)
        mis_z = float(constraint.get("mis_z") or 0)
        lof_z = float(constraint.get("lof_z") or 0)

        # Normalize: pLI is 0-1, oe ratios invert (lower = more constrained)
        # Z-scores can be negative; shift to 0-1 range for display
        metrics = {
            "pLI (LoF intolerance)": pli,
            "LoF constraint (1-LOEUF)": max(0, 1 - oe_lof),
            "Missense constraint (1-o/e)": max(0, 1 - oe_mis),
            "Missense Z": min(1, max(0, (mis_z + 2) / 6)),  # normalize -2..4 → 0..1
            "LoF Z": min(1, max(0, (lof_z + 2) / 6)),
        }

        categories = list(metrics.keys())
        values = list(metrics.values())
        # Close the polygon
        categories.append(categories[0])
        values.append(values[0])

        fig = go.Figure()
        fig.add_trace(go.Scatterpolar(
            r=values,
            theta=categories,
            fill="toself",
            fillcolor=f"rgba(99, 110, 250, 0.3)",
            line=dict(color=COLORS["muted"], width=2),
            marker=dict(size=6, color=COLORS["primary"]),
            hovertemplate="<b>%{theta}</b><br>Score: %{r:.3f}<extra></extra>",
        ))

        fig.update_layout(
            **_base_layout(
                title=dict(text=f"{gene} — Gene Constraint (gnomAD v4)", font=dict(size=16)),
                height=420,
                polar=dict(
                    bgcolor=COLORS["card_bg"],
                    radialaxis=dict(
                        visible=True, range=[0, 1],
                        gridcolor=COLORS["grid"],
                        tickfont=dict(size=10),
                    ),
                    angularaxis=dict(
                        gridcolor=COLORS["grid"],
                        tickfont=dict(size=11),
                    ),
                ),
                showlegend=False,
            )
        )

        # Add annotation with raw values
        raw_text = (
            f"pLI={pli:.3f}  |  LOEUF={oe_lof:.2f}  |  "
            f"mis_z={mis_z:.2f}  |  lof_z={lof_z:.2f}"
        )
        fig.add_annotation(
            text=raw_text, xref="paper", yref="paper",
            x=0.5, y=-0.08, showarrow=False,
            font=dict(size=10, color=COLORS["text"]),
        )
        return fig
    except Exception as e:
        logger.warning(f"chart_gene_constraint error: {e}")
        return None


def chart_source_coverage(profile_json: str) -> go.Figure | None:
    """Donut chart showing record coverage across NeuroPlex sources.

    Args:
        profile_json: JSON string from query_neuroplex_gene_profile.
    """
    try:
        profile = json.loads(profile_json)
        sources = profile.get("sources", {})
        gene = profile.get("gene", "")

        # Extract counts per source
        source_labels = {
            "opentargets": ("OpenTargets", "n_associations"),
            "gnomad": ("gnomAD", "n_variants"),
            "kegg": ("KEGG", "n_pathways"),
            "cbioportal": ("cBioPortal", "n_mutations"),
            "uniprot": ("UniProt", "n_entries"),
            "ncbi_gtr": ("NCBI GTR", "n_tests"),
        }

        labels = []
        values = []
        for key, (display_name, count_key) in source_labels.items():
            data = sources.get(key, {})
            count = data.get(count_key, 0) if "error" not in data else 0
            if count > 0:
                labels.append(display_name)
                values.append(count)

        if not values:
            return None

        fig = go.Figure()
        fig.add_trace(go.Pie(
            labels=labels,
            values=values,
            hole=0.55,
            marker=dict(colors=DISEASE_PALETTE[:len(labels)]),
            textinfo="label+value",
            textfont=dict(size=12),
            hovertemplate="<b>%{label}</b><br>Records: %{value}<br>%{percent}<extra></extra>",
        ))

        fig.update_layout(
            **_base_layout(
                title=dict(text=f"{gene} — NeuroPlex Source Coverage", font=dict(size=16)),
                height=380,
                showlegend=True,
                legend=dict(orientation="h", yanchor="bottom", y=-0.15, xanchor="center", x=0.5),
                annotations=[dict(
                    text=f"<b>{sum(values)}</b><br>records",
                    x=0.5, y=0.5, font=dict(size=18, color=COLORS["text"]),
                    showarrow=False,
                )],
            )
        )
        return fig
    except Exception as e:
        logger.warning(f"chart_source_coverage error: {e}")
        return None


def chart_panel_disease_heatmap(landscape_json: str) -> go.Figure | None:
    """Heatmap of gene-disease association scores across the NeuroPlex panel.

    Args:
        landscape_json: JSON string from query_neuroplex_disease_landscape.
    """
    try:
        data = json.loads(landscape_json)
        results = data.get("results", [])
        if not results:
            return None

        # Build gene x disease matrix
        genes = sorted(set(r.get("gene_symbol", "") for r in results))
        diseases = []
        seen_diseases = set()
        for r in results:
            d = r.get("disease", "")
            if d and d not in seen_diseases:
                diseases.append(d)
                seen_diseases.add(d)
        diseases = diseases[:12]  # limit for readability

        # Build matrix
        matrix = []
        for disease in diseases:
            row = []
            for gene in genes:
                score = next(
                    (float(r.get("score") or 0)
                     for r in results
                     if r.get("gene_symbol") == gene and r.get("disease") == disease),
                    0.0
                )
                row.append(score)
            matrix.append(row)

        fig = go.Figure()
        fig.add_trace(go.Heatmap(
            z=matrix,
            x=genes,
            y=[d[:35] for d in diseases],
            colorscale=[
                [0, COLORS["card_bg"]],
                [0.3, "#1a3a5c"],
                [0.6, "#2a6fa8"],
                [1.0, COLORS["primary"]],
            ],
            hovertemplate="<b>%{x}</b> × %{y}<br>Score: %{z:.4f}<extra></extra>",
            colorbar=dict(title="Score", thickness=15),
        ))

        fig.update_layout(
            **_base_layout(
                title=dict(text="NeuroPlex Gene-Disease Association Heatmap", font=dict(size=16)),
                height=max(400, len(diseases) * 35 + 120),
                xaxis=dict(tickangle=-45, tickfont=dict(size=11)),
                yaxis=dict(tickfont=dict(size=11)),
            )
        )
        return fig
    except Exception as e:
        logger.warning(f"chart_panel_disease_heatmap error: {e}")
        return None


def chart_constraint_comparison(variant_summary_json: str) -> go.Figure | None:
    """Scatter plot comparing pLI vs LOEUF across all panel genes.

    Highlights genes under strong constraint (high pLI, low LOEUF).

    Args:
        variant_summary_json: JSON string from query_neuroplex_variant_summary.
    """
    try:
        data = json.loads(variant_summary_json)
        constraints = data.get("gene_constraint", [])
        if not constraints:
            return None

        genes = [c.get("gene_symbol", "") for c in constraints]
        pli_vals = [float(c.get("pLI") or 0) for c in constraints]
        loeuf_vals = [float(c.get("oe_lof_upper") or 1) for c in constraints]
        mis_z_vals = [float(c.get("mis_z") or 0) for c in constraints]

        # Size by absolute mis_z
        sizes = [max(8, min(30, abs(mz) * 10)) for mz in mis_z_vals]

        # Color by pLI (higher = more constrained = red)
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=loeuf_vals,
            y=pli_vals,
            mode="markers+text",
            marker=dict(
                size=sizes,
                color=pli_vals,
                colorscale=[[0, COLORS["secondary"]], [1, COLORS["primary"]]],
                showscale=True,
                colorbar=dict(title="pLI", thickness=12),
                line=dict(width=1, color=COLORS["grid"]),
            ),
            text=genes,
            textposition="top center",
            textfont=dict(size=9, color=COLORS["text"]),
            hovertemplate=(
                "<b>%{text}</b><br>"
                "pLI: %{y:.3f}<br>"
                "LOEUF: %{x:.3f}<br>"
                "<extra></extra>"
            ),
        ))

        # Reference lines
        fig.add_hline(y=0.9, line_dash="dash", line_color=COLORS["primary"],
                      annotation_text="pLI=0.9 (LoF intolerant)", annotation_position="top left")
        fig.add_vline(x=0.35, line_dash="dash", line_color=COLORS["warning"],
                      annotation_text="LOEUF=0.35", annotation_position="top right")

        fig.update_layout(
            **_base_layout(
                title=dict(text="NeuroPlex Panel — Gene Constraint (pLI vs LOEUF)", font=dict(size=16)),
                xaxis_title="LOEUF (LoF observed/expected upper bound)",
                yaxis_title="pLI (probability of LoF intolerance)",
                height=500,
                showlegend=False,
            )
        )
        return fig
    except Exception as e:
        logger.warning(f"chart_constraint_comparison error: {e}")
        return None

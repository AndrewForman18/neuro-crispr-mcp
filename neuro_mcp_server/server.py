"""Neuro CRISPR MCP Server — CRISPR Atlas + CRISPRbrain tools.

Provides tools to explore:
- PerturbAI Whole Brain CRISPR Atlas (MOUSE, in vivo brain, 7.7M cells, 2,046 KOs)
- CRISPRbrain (HUMAN, iPSC-derived neurons/microglia/astrocytes, 127 screens)
- Human expression (BioFINDER, ROSMAP cohorts)
"""

import json
import logging
import os
from typing import Any

from databricks.sdk import WorkspaceClient

from config.neuroplex_config import load_config

logger = logging.getLogger(__name__)

# ─── Configuration ───────────────────────────────────────────────────────────
CFG = load_config()
CATALOG = CFG.catalog
SCHEMA = CFG.schema
# Mouse (CRISPR Atlas)
TABLE_ATLAS = f"{CATALOG}.{SCHEMA}.wholebrain_crispr_atlas"
TABLE_DIFF_EXPR = f"{CATALOG}.{SCHEMA}.crispr_atlas_diff_expr"
TABLE_CELL_META = f"{CATALOG}.{SCHEMA}.crispr_atlas_cell_metadata"
TABLE_GENE_META = f"{CATALOG}.{SCHEMA}.crispr_atlas_gene_metadata"
TABLE_NDD_GENES = f"{CATALOG}.{SCHEMA}.crispr_atlas_ndd_genes"
# Human (CRISPRbrain)
TABLE_CRISPRBRAIN = f"{CATALOG}.{SCHEMA}.crisprbrain_screens"
# Human expression
TABLE_HUMAN_EXPR = f"{CATALOG}.{SCHEMA}.gene_expression_matrix"
# GTEx brain expression (baseline)
TABLE_GTEX = f"{CATALOG}.{SCHEMA}.gtex_brain_expression"
# LINCS L1000 drug perturbation signatures
TABLE_LINCS = f"{CATALOG}.{SCHEMA}.lincs_l1000_signatures"
# ChEMBL orexin receptor pharmacology
TABLE_CHEMBL = f"{CATALOG}.{SCHEMA}.chembl_orexin_pharmacology"
# Tahoe-100M drug perturbation atlas
TABLE_TAHOE = f"{CATALOG}.{SCHEMA}.tahoe_100m"
TABLE_TAHOE_CLUSTERED = f"{CATALOG}.{SCHEMA}.tahoe_100m_clustered"
TABLE_TAHOE_GENE_VOCAB = f"{CATALOG}.{SCHEMA}.tahoe_100m_gene_vocab"

QUERY_TIMEOUT_SECONDS = 90

# SQL warehouse
_wh_id = os.environ.get("DATABRICKS_SQL_WAREHOUSE_HTTP_PATH", CFG.sql_warehouse_id)
SQL_WAREHOUSE_ID = _wh_id.split("/")[-1] if "/" in _wh_id else _wh_id

_ws_client = WorkspaceClient()


# ─── Database helper ─────────────────────────────────────────────────────────
def _execute_query(query: str, params: dict | None = None, timeout: int | None = None) -> list[dict[str, Any]]:
    """Execute SQL via Statement Execution API."""
    effective_timeout = timeout or QUERY_TIMEOUT_SECONDS
    logger.info(f"SQL ({len(query)} chars), timeout={effective_timeout}s")

    final_query = query
    if params:
        for key, value in params.items():
            placeholder = f":{key}"
            if isinstance(value, str):
                escaped = value.replace("'", "''")
                final_query = final_query.replace(placeholder, f"'{escaped}'")
            else:
                final_query = final_query.replace(placeholder, str(value))

    try:
        result = _ws_client.statement_execution.execute_statement(
            warehouse_id=SQL_WAREHOUSE_ID,
            statement=final_query,
            wait_timeout=f"{min(effective_timeout, 50)}s",
        )
        if result.status.state.value == "FAILED":
            raise RuntimeError(f"SQL failed: {result.status.error}")
        if not result.manifest or not result.result:
            return []
        columns = [col.name for col in result.manifest.schema.columns]
        rows = result.result.data_array or []
        logger.info(f"Query returned {len(rows)} rows")
        return [dict(zip(columns, row)) for row in rows]
    except Exception as e:
        logger.error(f"SQL FAILED: {e}")
        raise


# ─── MOUSE Tools (CRISPR Atlas) ─────────────────────────────────────────────

def query_knockout_effects(gene_target: str, cell_type: str = "", pval_threshold: float = 0.05, limit: int = 50) -> str:
    """Find genes affected when a specific gene is knocked out in MOUSE brain (CRISPR Atlas).

    Args:
        gene_target: Mouse gene knocked out (Title Case, e.g. Psen1, App, Hcrt).
        cell_type: Optional cell type filter (e.g. Glut, GABA, CTX).
        pval_threshold: Adjusted p-value threshold (default 0.05).
        limit: Max results (default 50, max 200).
    """
    limit = min(max(1, limit), 200)
    where = ["gene_target = :gene_target", "pvalue_adj < :pval"]
    params: dict[str, Any] = {"gene_target": gene_target, "pval": pval_threshold}
    if cell_type:
        where.append("group_name LIKE '%' || :cell_type || '%'")
        params["cell_type"] = cell_type

    query = f"""
        SELECT de_gene AS affected_gene, group_name AS cell_type,
               log2fc, pvalue_adj,
               n_cells_ko AS n_perturbed, n_cells_ctrl AS n_control
        FROM {TABLE_DIFF_EXPR}
        WHERE {' AND '.join(where)}
        ORDER BY ABS(log2fc) DESC
        LIMIT {limit}
    """
    results = _execute_query(query, params)
    return json.dumps({"species": "mouse", "dataset": "CRISPR Atlas (in vivo brain)", "knockout": gene_target, "cell_type_filter": cell_type or "all", "significant_effects": len(results), "results": results}, default=str)


def find_knockouts_affecting_gene(gene_name: str, cell_type: str = "", pval_threshold: float = 0.05, limit: int = 50) -> str:
    """Find which CRISPR knockouts affect expression of a gene in MOUSE brain.

    Args:
        gene_name: Downstream mouse gene (Title Case, e.g. Hcrt, Grin1, Drd2).
        cell_type: Optional cell type filter.
        pval_threshold: FDR threshold.
        limit: Max results.
    """
    limit = min(max(1, limit), 200)
    where = ["de_gene = :gene_name", "pvalue_adj < :pval"]
    params: dict[str, Any] = {"gene_name": gene_name, "pval": pval_threshold}
    if cell_type:
        where.append("group_name LIKE '%' || :cell_type || '%'")
        params["cell_type"] = cell_type

    query = f"""
        SELECT gene_target AS knockout, group_name AS cell_type,
               log2fc, pvalue_adj,
               n_cells_ko AS n_perturbed, n_cells_ctrl AS n_control
        FROM {TABLE_DIFF_EXPR}
        WHERE {' AND '.join(where)}
        ORDER BY ABS(log2fc) DESC
        LIMIT {limit}
    """
    results = _execute_query(query, params)
    return json.dumps({"species": "mouse", "dataset": "CRISPR Atlas", "target_gene": gene_name, "knockouts_affecting_gene": len(results), "results": results}, default=str)


def list_ndd_risk_genes(disorder: str = "NDD", fdr_threshold: float = 0.1, limit: int = 50) -> str:
    """List genes associated with neurodevelopmental disorders (ASD, DD, NDD).

    Args:
        disorder: ASD, DD, or NDD (default NDD).
        fdr_threshold: FDR threshold.
        limit: Max results.
    """
    limit = min(max(1, limit), 200)
    disorder = disorder.upper()
    fdr_col = f"FDR_TADA_{disorder}" if disorder in ("ASD", "DD", "NDD") else None
    if not fdr_col:
        return json.dumps({"error": f"Unsupported disorder. Use ASD, DD, or NDD."})

    query = f"""
        SELECT DISTINCT gene, gene_id, chromosome,
               {fdr_col} AS fdr, p_TADA_{disorder} AS pvalue,
               ASD72, DD309, NDD373, SCZ244
        FROM {TABLE_NDD_GENES}
        WHERE {fdr_col} IS NOT NULL AND {fdr_col} < :fdr_threshold
        ORDER BY {fdr_col} ASC
        LIMIT {limit}
    """
    results = _execute_query(query, {"fdr_threshold": fdr_threshold})
    return json.dumps({"disorder": disorder, "fdr_threshold": fdr_threshold, "risk_genes_found": len(results), "results": results}, default=str)


def list_brain_cell_types(region_filter: str = "", limit: int = 50) -> str:
    """List brain cell types in the MOUSE CRISPR Atlas."""
    limit = min(max(1, limit), 100)
    where = "WHERE passes_qc = true"
    params: dict[str, Any] = {}
    if region_filter:
        where += " AND (LOWER(predicted_class) LIKE LOWER('%' || :filt || '%') OR LOWER(region_level1) LIKE LOWER('%' || :filt || '%') OR LOWER(neuron_type) LIKE LOWER('%' || :filt || '%'))"
        params["filt"] = region_filter

    query = f"""
        SELECT predicted_class, predicted_subclass, neuron_type,
               region_level1, region_level2, COUNT(*) AS n_cells
        FROM {TABLE_CELL_META}
        {where}
        GROUP BY predicted_class, predicted_subclass, neuron_type, region_level1, region_level2
        ORDER BY n_cells DESC
        LIMIT {limit}
    """
    results = _execute_query(query, params if params else None, timeout=60)
    return json.dumps({"species": "mouse", "filter": region_filter or "all", "cell_types_found": len(results), "results": results}, default=str)


def list_knockout_targets(search: str = "", limit: int = 50) -> str:
    """List available CRISPR knockout targets in MOUSE atlas."""
    limit = min(max(1, limit), 200)
    params: dict[str, Any] = {}
    where = "WHERE gene_target IS NOT NULL AND gene_target != 'non-targeting'"
    if search:
        where += " AND LOWER(gene_target) LIKE LOWER('%' || :search || '%')"
        params["search"] = search

    query = f"""
        SELECT DISTINCT gene_target
        FROM {TABLE_CELL_META}
        {where}
        ORDER BY gene_target
        LIMIT {limit}
    """
    results = _execute_query(query, params if params else None, timeout=30)
    return json.dumps({"species": "mouse", "search": search or "all", "targets_found": len(results), "results": results}, default=str)


# ─── HUMAN Tools (CRISPRbrain) ──────────────────────────────────────────────

def query_human_screen_hits(gene: str, cell_type: str = "", phenotype: str = "", limit: int = 50) -> str:
    """Query HUMAN CRISPRbrain screen results for a gene.

    Returns phenotype scores, p-values, and hit classification from iPSC-derived
    human neurons, microglia, and astrocytes.

    Args:
        gene: Human gene symbol (UPPER CASE, e.g. PSEN1, APP, MAPT, SOD2).
        cell_type: Optional cell type filter (e.g. Glutamatergic, Microglia, Astrocyte).
        phenotype: Optional phenotype filter (e.g. Survival, Tau, TDP-43).
        limit: Max results (default 50).
    """
    limit = min(max(1, limit), 200)
    where = ["UPPER(gene) = UPPER(:gene)"]
    params: dict[str, Any] = {"gene": gene}
    if cell_type:
        where.append("LOWER(cell_type) LIKE LOWER('%' || :cell_type || '%')")
        params["cell_type"] = cell_type
    if phenotype:
        where.append("LOWER(phenotype_name) LIKE LOWER('%' || :phenotype || '%')")
        params["phenotype"] = phenotype

    query = f"""
        SELECT screen_name, gene, cell_type, genotype, crispr_mode,
               phenotype_name, phenotype_score, pvalue, gene_score, hit_class
        FROM {TABLE_CRISPRBRAIN}
        WHERE {' AND '.join(where)}
        ORDER BY ABS(phenotype_score) DESC
        LIMIT {limit}
    """
    results = _execute_query(query, params)
    return json.dumps({"species": "human", "dataset": "CRISPRbrain (iPSC-derived cells)", "gene": gene, "screens_found": len(results), "results": results}, default=str)


def list_human_screens(cell_type: str = "", limit: int = 50) -> str:
    """List available HUMAN CRISPRbrain screens with hit counts.

    Args:
        cell_type: Optional filter (e.g. Glutamatergic, Microglia, Astrocyte).
        limit: Max results.
    """
    limit = min(max(1, limit), 100)
    where = ""
    params: dict[str, Any] = {}
    if cell_type:
        where = "WHERE LOWER(cell_type) LIKE LOWER('%' || :cell_type || '%')"
        params["cell_type"] = cell_type

    query = f"""
        SELECT screen_name, cell_type, genotype, crispr_mode, phenotype_name,
               COUNT(*) AS n_genes,
               SUM(CASE WHEN hit_class != 'none' THEN 1 ELSE 0 END) AS n_hits
        FROM {TABLE_CRISPRBRAIN}
        {where}
        GROUP BY screen_name, cell_type, genotype, crispr_mode, phenotype_name
        ORDER BY n_hits DESC
        LIMIT {limit}
    """
    results = _execute_query(query, params if params else None)
    return json.dumps({"species": "human", "dataset": "CRISPRbrain", "screens_found": len(results), "results": results}, default=str)


def compare_cross_species(gene: str, limit: int = 30) -> str:
    """Compare a gene's CRISPR effects across MOUSE (Atlas) and HUMAN (CRISPRbrain).

    Queries both datasets to show cross-species conservation of gene function.

    Args:
        gene: Gene symbol (auto-converted: UPPER for human, Title for mouse).
        limit: Max results per species.
    """
    limit = min(max(1, limit), 100)
    human_gene = gene.upper()
    mouse_gene = gene.capitalize()

    # Human CRISPRbrain
    human_query = f"""
        SELECT screen_name, cell_type, phenotype_name, phenotype_score, pvalue, hit_class
        FROM {TABLE_CRISPRBRAIN}
        WHERE UPPER(gene) = '{human_gene}'
        ORDER BY ABS(phenotype_score) DESC
        LIMIT {limit}
    """
    human_results = _execute_query(human_query, timeout=30)

    # Mouse CRISPR Atlas (as KO target)
    mouse_query = f"""
        SELECT group_name AS cell_type, COUNT(*) AS n_sig_effects,
               ROUND(AVG(ABS(log2fc)), 3) AS avg_abs_log2fc,
               MIN(pvalue_adj) AS best_pval
        FROM {TABLE_DIFF_EXPR}
        WHERE gene_target = '{mouse_gene}' AND pvalue_adj < 0.05
        GROUP BY group_name
        ORDER BY n_sig_effects DESC
        LIMIT {limit}
    """
    mouse_results = _execute_query(mouse_query, timeout=30)

    return json.dumps({
        "gene": gene,
        "human": {"dataset": "CRISPRbrain (iPSC neurons/glia)", "gene_symbol": human_gene, "screens_with_data": len(human_results), "results": human_results},
        "mouse": {"dataset": "CRISPR Atlas (in vivo brain)", "gene_symbol": mouse_gene, "cell_types_affected": len(mouse_results), "results": mouse_results},
    }, default=str)


# ─── Pharmacology Tools ──────────────────────────────────────────────────────

def query_baseline_expression(gene_symbol: str, tissue: str = "") -> str:
    """Query baseline gene expression across brain regions from GTEx.

    Returns median TPM expression levels per brain subregion.
    Useful for understanding where orexin system genes (HCRT, HCRTR1, HCRTR2)
    are normally expressed in the brain.

    Args:
        gene_symbol: Gene symbol (UPPER CASE, e.g. HCRT, HCRTR1, HCRTR2, APP).
        tissue: Optional brain region filter (e.g. Hypothalamus, Cortex, Cerebellum).
    """
    where = ["UPPER(gene_symbol) = UPPER(:gene_symbol)"]
    params: dict[str, Any] = {"gene_symbol": gene_symbol}
    if tissue:
        where.append("LOWER(tissue) LIKE LOWER('%' || :tissue || '%')")
        params["tissue"] = tissue

    query = f"""
        SELECT gene_symbol, tissue, tissue_id, median_tpm, n_samples, source
        FROM {TABLE_GTEX}
        WHERE {' AND '.join(where)}
        ORDER BY median_tpm DESC
    """
    results = _execute_query(query, params, timeout=30)
    if not results:
        return json.dumps({"error": f"Gene '{gene_symbol}' not found in GTEx brain expression data."})
    return json.dumps({
        "species": "human",
        "dataset": "GTEx v8 (brain regions)",
        "gene": gene_symbol,
        "brain_regions": len(results),
        "results": results,
    }, default=str)


def query_drug_perturbations(compound: str = "", gene_symbol: str = "", cell_line: str = "", neuronal_only: bool = False, limit: int = 50) -> str:
    """Query LINCS L1000 drug perturbation signatures for orexin-related compounds.

    Returns gene expression z-scores showing how compounds affect transcription
    in neuronal cell lines.

    Args:
        compound: Drug name (e.g. suvorexant, lemborexant, daridorexant, SB-334867).
        gene_symbol: Optional gene to check specific effect on (e.g. HCRT, GABRA1).
        cell_line: Optional cell line filter (e.g. SHSY5Y, NPC, NEU).
        neuronal_only: If True, restrict to neuronal lines only.
        limit: Max results (default 50).
    """
    limit = min(max(1, limit), 200)
    where: list[str] = []
    params: dict[str, Any] = {}

    if compound:
        where.append("LOWER(compound) = LOWER(:compound)")
        params["compound"] = compound
    if gene_symbol:
        where.append("UPPER(gene_symbol) = UPPER(:gene_symbol)")
        params["gene_symbol"] = gene_symbol
    if cell_line:
        where.append("UPPER(cell_line) = UPPER(:cell_line)")
        params["cell_line"] = cell_line
    if neuronal_only:
        where.append("is_neuronal = true")

    if not where:
        return json.dumps({"error": "Provide at least one of: compound, gene_symbol, or cell_line."})

    query = f"""
        SELECT compound, moa, cell_line, is_neuronal, gene_symbol,
               zscore, dose_um, timepoint_h, approval_status, brand_name
        FROM {TABLE_LINCS}
        WHERE {' AND '.join(where)}
        ORDER BY ABS(zscore) DESC
        LIMIT {limit}
    """
    results = _execute_query(query, params, timeout=30)
    return json.dumps({
        "dataset": "LINCS L1000 (drug perturbation signatures)",
        "compound_filter": compound or "all",
        "gene_filter": gene_symbol or "all",
        "cell_line_filter": cell_line or ("neuronal only" if neuronal_only else "all"),
        "signatures_found": len(results),
        "results": results,
    }, default=str)


def query_receptor_pharmacology(compound: str = "", target: str = "", assay_type: str = "", limit: int = 50) -> str:
    """Query ChEMBL binding/activity data for orexin receptor ligands.

    Returns IC50, Ki, Kd, and pChEMBL values for compounds tested against
    OX1R (HCRTR1) and OX2R (HCRTR2).

    Args:
        compound: Drug name (e.g. suvorexant, lemborexant, daridorexant).
        target: Receptor filter (e.g. OX1R, OX2R, HCRTR1, HCRTR2).
        assay_type: Activity type filter (e.g. Ki, IC50, Kd).
        limit: Max results.
    """
    limit = min(max(1, limit), 200)
    where: list[str] = []
    params: dict[str, Any] = {}

    if compound:
        where.append("LOWER(molecule_name) LIKE LOWER('%' || :compound || '%')")
        params["compound"] = compound
    if target:
        target_upper = target.upper()
        where.append("(UPPER(target_alias) = :target OR UPPER(target_gene) = :target)")
        params["target"] = target_upper
    if assay_type:
        where.append("UPPER(standard_type) = UPPER(:assay_type)")
        params["assay_type"] = assay_type

    if not where:
        # Default: show all clinical compounds
        where.append("molecule_name IS NOT NULL AND molecule_name != ''")

    query = f"""
        SELECT molecule_name, molecule_chembl_id, target_gene, target_alias,
               compound_brand, compound_company, selectivity_class, year_approved,
               standard_type, standard_value, standard_units, pchembl_value,
               assay_description
        FROM {TABLE_CHEMBL}
        WHERE {' AND '.join(where)}
        ORDER BY pchembl_value DESC NULLS LAST
        LIMIT {limit}
    """
    results = _execute_query(query, params, timeout=30)
    return json.dumps({
        "dataset": "ChEMBL (orexin receptor pharmacology)",
        "compound_filter": compound or "all",
        "target_filter": target or "OX1R + OX2R",
        "activities_found": len(results),
        "results": results,
    }, default=str)


# ─── Shared Tools ────────────────────────────────────────────────────────────

def query_human_expression(gene_symbol: str, cohort: str = "", tissue: str = "") -> str:
    """Query HUMAN gene expression from BioFINDER/ROSMAP neurology cohorts."""
    params: dict[str, Any] = {"gene_symbol": gene_symbol}
    where = ["LOWER(gene_symbol) = LOWER(:gene_symbol)"]
    if cohort:
        where.append("LOWER(cohort) = LOWER(:cohort)")
        params["cohort"] = cohort
    if tissue:
        where.append("LOWER(tissue) LIKE LOWER('%' || :tissue || '%')")
        params["tissue"] = tissue

    query = f"""
        SELECT gene_symbol, cohort, tissue,
               COUNT(*) AS n_samples,
               ROUND(AVG(tpm), 4) AS mean_tpm,
               ROUND(PERCENTILE(tpm, 0.5), 4) AS median_tpm,
               ROUND(STDDEV(tpm), 4) AS std_tpm
        FROM {TABLE_HUMAN_EXPR}
        WHERE {' AND '.join(where)}
        GROUP BY gene_symbol, cohort, tissue
        ORDER BY mean_tpm DESC
    """
    results = _execute_query(query, params, timeout=60)
    if not results:
        return json.dumps({"error": f"Gene '{gene_symbol}' not found in human expression matrix."})
    return json.dumps({"species": "human", "dataset": "gene_expression_matrix (BioFINDER, ROSMAP)", "gene": gene_symbol, "expression_summary": results}, default=str)


def summarize_atlas() -> str:
    """Get high-level summary of BOTH datasets (Mouse Atlas + Human CRISPRbrain)."""
    mouse_targets = _execute_query(f"SELECT COUNT(DISTINCT gene_target) AS n FROM {TABLE_CELL_META} WHERE gene_target IS NOT NULL AND gene_target != 'non-targeting'", timeout=30)
    human_screens = _execute_query(f"SELECT COUNT(DISTINCT screen_name) AS n, COUNT(*) AS total_rows FROM {TABLE_CRISPRBRAIN}", timeout=30)

    return json.dumps({
        "mouse_crispr_atlas": {
            "source": "PerturbAI Whole Brain CRISPR Atlas",
            "species": "mouse (in vivo brain)",
            "knockout_targets": mouse_targets[0]["n"] if mouse_targets else "unknown",
            "diff_expr_table": TABLE_DIFF_EXPR,
            "description": "745M pre-computed DE results across hundreds of brain cell types",
        },
        "human_crisprbrain": {
            "source": "CRISPRbrain.org (Kampmann Lab, UCSF)",
            "species": "human (iPSC-derived neurons, microglia, astrocytes)",
            "screens": human_screens[0]["n"] if human_screens else "unknown",
            "total_rows": human_screens[0]["total_rows"] if human_screens else "unknown",
            "phenotypes": "Survival, ROS, iron, lipids, lysosomes, tau aggregation, TDP-43, phagocytosis",
            "table": TABLE_CRISPRBRAIN,
        },
        "human_expression": {
            "source": "BioFINDER, ROSMAP",
            "species": "human",
            "data": "Full transcriptome RNA-seq (brain + blood)",
            "table": TABLE_HUMAN_EXPR,
        },
    }, default=str)


# ─── Tahoe-100M Tools ────────────────────────────────────────────

def list_tahoe_drugs(search: str = "", cell_line: str = "", limit: int = 50) -> str:
    """List compounds in the Tahoe-100M drug perturbation atlas with MOA and cell counts.

    Args:
        search: Optional drug name substring filter.
        cell_line: Optional cell line filter.
        limit: Max results (default 50).
    """
    limit = min(max(1, limit), 200)
    where_clauses = ["drug IS NOT NULL"]
    params: dict[str, Any] = {}
    if search:
        where_clauses.append("LOWER(drug) LIKE LOWER('%' || :search || '%')")
        params["search"] = search
    if cell_line:
        where_clauses.append("LOWER(cell_line_id) LIKE LOWER('%' || :cell_line || '%')")
        params["cell_line"] = cell_line
    where = "WHERE " + " AND ".join(where_clauses)
    query = f"""
        SELECT drug, moa,
               COUNT(*) AS n_cells,
               COUNT(DISTINCT cell_line_id) AS n_cell_lines
        FROM {TABLE_TAHOE_CLUSTERED}
        {where}
        GROUP BY drug, moa
        ORDER BY n_cells DESC
        LIMIT {limit}
    """
    results = _execute_query(query, params if params else None, timeout=60)
    return json.dumps({
        "dataset": "Tahoe-100M Drug Perturbation Atlas",
        "search_filter": search or "all",
        "cell_line_filter": cell_line or "all",
        "drugs_found": len(results),
        "results": results,
    }, default=str)


def query_tahoe_drug_clusters(drug: str = "", cell_line: str = "", limit: int = 50) -> str:
    """Query cluster distribution for a drug in Tahoe-100M.

    Shows how cells treated with a compound distribute across transcriptomic clusters
    and cell lines.

    Args:
        drug: Drug name (or substring) to query.
        cell_line: Optional cell line filter.
        limit: Max results.
    """
    limit = min(max(1, limit), 200)
    where_clauses: list[str] = []
    params: dict[str, Any] = {}
    if drug:
        where_clauses.append("LOWER(drug) LIKE LOWER('%' || :drug || '%')")
        params["drug"] = drug
    if cell_line:
        where_clauses.append("LOWER(cell_line_id) LIKE LOWER('%' || :cell_line || '%')")
        params["cell_line"] = cell_line
    where = ("WHERE " + " AND ".join(where_clauses)) if where_clauses else ""
    query = f"""
        SELECT drug, cell_line_id, moa, cluster,
               COUNT(*) AS n_cells
        FROM {TABLE_TAHOE_CLUSTERED}
        {where}
        GROUP BY drug, cell_line_id, moa, cluster
        ORDER BY drug, n_cells DESC
        LIMIT {limit}
    """
    results = _execute_query(query, params if params else None, timeout=60)
    return json.dumps({
        "dataset": "Tahoe-100M (clustered)",
        "drug_filter": drug or "all",
        "cell_line_filter": cell_line or "all",
        "records": len(results),
        "results": results,
    }, default=str)


def query_tahoe_moa_distribution(moa_filter: str = "", limit: int = 50) -> str:
    """Summarize Tahoe-100M compounds by mechanism of action (MOA).

    Shows how many drugs and cells belong to each MOA class.

    Args:
        moa_filter: Optional MOA substring filter.
        limit: Max results.
    """
    limit = min(max(1, limit), 200)
    params: dict[str, Any] = {}
    where = "WHERE moa IS NOT NULL"
    if moa_filter:
        where += " AND LOWER(moa) LIKE LOWER('%' || :moa_filter || '%')"
        params["moa_filter"] = moa_filter
    query = f"""
        SELECT moa,
               COUNT(DISTINCT drug) AS n_drugs,
               COUNT(*) AS n_cells,
               COUNT(DISTINCT cell_line_id) AS n_cell_lines
        FROM {TABLE_TAHOE_CLUSTERED}
        {where}
        GROUP BY moa
        ORDER BY n_drugs DESC
        LIMIT {limit}
    """
    results = _execute_query(query, params if params else None, timeout=60)
    return json.dumps({
        "dataset": "Tahoe-100M",
        "moa_filter": moa_filter or "all",
        "moa_classes_found": len(results),
        "results": results,
    }, default=str)


# ─── NeuroPlex Tools (Multi-Source Gene Intelligence) ───────────────────────

from neuro_mcp_server.neuroplex_query import (
    query_neuroplex_gene_profile as _npx_gene_profile,
    query_neuroplex_disease_landscape as _npx_disease_landscape,
    query_neuroplex_variant_summary as _npx_variant_summary,
    chart_disease_associations,
    chart_variant_landscape,
    chart_gene_constraint,
    chart_source_coverage,
    chart_panel_disease_heatmap,
    chart_constraint_comparison,
)


def query_neuroplex_gene(gene_symbol: str) -> str:
    """Query unified NeuroPlex gene profile across 6 biomedical sources.

    Aggregates intelligence from OpenTargets (disease associations), gnomAD
    (ClinVar variants + gene constraint), KEGG (pathways), cBioPortal
    (somatic mutations), UniProt (protein annotations), and NCBI GTR
    (genetic tests) for a single gene in the neuroscience target panel.

    Args:
        gene_symbol: Gene symbol (UPPER CASE, e.g. HCRT, PSEN1, APP, SOD1, MAPT).
    """
    return _npx_gene_profile(gene_symbol, _execute_query_fn=_execute_query)


def query_neuroplex_diseases(disease_filter: str = "", limit: int = 15) -> str:
    """Query NeuroPlex disease-gene landscape from OpenTargets.

    Returns gene-disease association scores, optionally filtered by disease name.
    Shows which panel genes are linked to a specific neurological condition.

    Args:
        disease_filter: Disease name substring (e.g. narcolepsy, Alzheimer, Parkinson).
        limit: Max diseases per gene (default 15).
    """
    return _npx_disease_landscape(disease_filter, limit, _execute_query_fn=_execute_query)


def query_neuroplex_variants(gene_symbol: str = "") -> str:
    """Query NeuroPlex variant summary from gnomAD v4.

    Returns clinical significance distribution and gene constraint metrics
    (pLI, LOEUF, missense Z-score) across the neuroscience target panel.

    Args:
        gene_symbol: Optional gene filter (UPPER CASE). If empty, returns panel-wide summary.
    """
    return _npx_variant_summary(gene_symbol, _execute_query_fn=_execute_query)


def get_protein_structure(gene_symbol: str, uniprot_accession: str = "") -> str:
    """Get AlphaFold2 predicted protein structure from the EBI AlphaFold database.

    Queries the EBI AlphaFold protein structure database for a human protein.
    Returns structure metadata and direct PDB/CIF download URLs. Per-residue
    pLDDT confidence scores are embedded in the PDB file:
    >90 = very high, 70-90 = confident, 50-70 = low, <50 = disordered/flexible.

    Args:
        gene_symbol: Human gene symbol (UPPER CASE, e.g. PSEN1, APP, MAPT, HCRT, LRRK2, SOD1).
        uniprot_accession: Optional UniProt accession (e.g. P49768 for PSEN1). Skips
                           the UniProt lookup step when provided.
    """
    import urllib.request
    import urllib.parse

    try:
        accession = uniprot_accession.strip().upper() if uniprot_accession else ""
        resolved_gene = gene_symbol.upper()
        protein_name = ""
        seq_len = None

        # Step 1: resolve UniProt accession if not supplied
        if not accession:
            for query in [
                f"gene_exact:{urllib.parse.quote(gene_symbol)}+AND+organism_id:9606+AND+reviewed:true",
                f"gene:{urllib.parse.quote(gene_symbol)}+AND+organism_id:9606+AND+reviewed:true",
            ]:
                url = (
                    f"https://rest.uniprot.org/uniprotkb/search?query={query}"
                    f"&fields=accession,gene_names,protein_name,sequence&format=json&size=1"
                )
                with urllib.request.urlopen(url, timeout=10) as resp:
                    data = json.loads(resp.read())
                hits = data.get("results", [])
                if hits:
                    entry = hits[0]
                    accession = entry["primaryAccession"]
                    gene_list = entry.get("genes", [])
                    if gene_list:
                        resolved_gene = gene_list[0].get("geneName", {}).get("value", gene_symbol.upper())
                    pn_block = entry.get("proteinDescription", {}).get("recommendedName", {})
                    protein_name = pn_block.get("fullName", {}).get("value", "")
                    seq_len = entry.get("sequence", {}).get("length")
                    break

            if not accession:
                return json.dumps({
                    "error": (
                        f"No reviewed Swiss-Prot entry found for '{gene_symbol}' (human, organism_id 9606). "
                        f"Try supplying uniprot_accession directly to bypass this lookup."
                    ),
                    "gene_symbol": gene_symbol,
                })

        # Step 2: query EBI AlphaFold REST API
        ebi_url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
        with urllib.request.urlopen(ebi_url, timeout=15) as resp:
            predictions = json.loads(resp.read())

        if not predictions:
            return json.dumps({
                "error": (
                    f"No AlphaFold structure found for accession '{accession}'. "
                    f"The protein may not be covered by the AlphaFold database."
                ),
                "uniprot_accession": accession,
                "gene_symbol": resolved_gene,
            })

        pred = predictions[0]
        return json.dumps({
            "source": "EBI AlphaFold Database",
            "gene_symbol": resolved_gene,
            "uniprot_accession": accession,
            "protein_name": protein_name or pred.get("uniprotDescription", ""),
            "sequence_length": seq_len or pred.get("seqRes"),
            "alphafold_model_version": pred.get("latestVersion"),
            "model_created_date": pred.get("modelCreatedDate"),
            "last_updated": pred.get("lastUpdated"),
            "pdb_url": pred.get("pdbUrl", ""),
            "cif_url": pred.get("cifUrl", pred.get("bcifUrl", "")),
            "pae_image_url": pred.get("paeImageUrl", ""),
            "pae_doc_url": pred.get("paeDocUrl", ""),
            "alphafold_page_url": f"https://alphafold.ebi.ac.uk/entry/{accession}",
            "plddt_note": (
                "Per-residue pLDDT confidence scores are embedded in the PDB file. "
                ">90=very high confidence, 70-90=confident, 50-70=low, <50=disordered/flexible."
            ),
        }, default=str)

    except Exception as exc:
        logger.error(f"get_protein_structure failed for gene={gene_symbol!r}, accession={uniprot_accession!r}: {exc}")
        return json.dumps({
            "error": str(exc),
            "gene_symbol": gene_symbol,
            "uniprot_accession": uniprot_accession,
        })

"""Neuro CRISPR MCP — Streamlit Chat UI with Species Toggle.

A ReAct agent exploring CRISPR perturbation data across species:
- MOUSE: PerturbAI Whole Brain CRISPR Atlas (7.7M cells, 2,046 KOs)
- HUMAN: CRISPRbrain (iPSC neurons/microglia/astrocytes, 127 screens)
- HUMAN: BioFINDER/ROSMAP expression (full transcriptome, brain + blood)
"""

import json
import logging
import os
import sys
import time
import traceback
import uuid
from datetime import datetime

import streamlit as st

from config.neuroplex_config import load_config

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s", stream=sys.stdout)
logger = logging.getLogger("neuro_crispr_mcp")
CFG = load_config()

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from neuro_mcp_server.server import (
    _execute_query,
    query_knockout_effects,
    find_knockouts_affecting_gene,
    list_ndd_risk_genes,
    list_brain_cell_types,
    list_knockout_targets,
    query_human_screen_hits,
    list_human_screens,
    compare_cross_species,
    query_human_expression,
    query_baseline_expression,
    query_drug_perturbations,
    query_receptor_pharmacology,
    summarize_atlas,
    list_tahoe_drugs,
    query_tahoe_drug_clusters,
    query_tahoe_moa_distribution,
    query_neuroplex_gene,
    query_neuroplex_diseases,
    query_neuroplex_variants,
)
from neuro_mcp_server.neuroplex_query import (
    chart_disease_associations,
    chart_variant_landscape,
    chart_gene_constraint,
    chart_source_coverage,
    chart_panel_disease_heatmap,
    chart_constraint_comparison,
)

# ── i18n Translations ──
_I18N = {
    "en": {
        "subtitle": "Cross-species perturbation explorer",
        "sql_ready": "SQL Warehouse: {msg}",
        "sql_error": "SQL: {msg}",
        "dataset_focus": "\U0001F50D Dataset Focus",
        "btn_all": "All",
        "btn_none": "None",
        "ds_mouse": "\U0001F42D Mouse (CRISPR Atlas + IMPC)",
        "ds_human": "\U0001F9EC Human CRISPRbrain",
        "ds_crispr": "\U0001F9EC CRISPR Perturbation (Mouse + Human)",
        "ds_short_crispr": "CRISPR",
        "ds_expression": "\U0001F4CA Expression (GTEx + Cohorts)",
        "ds_pharma": "\U0001F48A Pharmacology (LINCS + ChEMBL)",
        "ds_tahoe": "\U0001F52C Tahoe-100M Drug Atlas",
        "ds_genetics": "\U0001F9EC Genetics (gnomAD + CIViC + GTR)",
        "ds_druggability": "\U0001F3AF Druggability (OpenTargets)",
        "ds_protein": "\U0001F52C Protein (HPA + Monarch)",
        "ds_pathways": "\U0001F5FA\uFE0F Pathways (UniProt + Reactome + KEGG)",
        "ds_short_genetics": "Genetics",
        "ds_short_druggability": "Druggability",
        "ds_short_protein": "Protein",
        "ds_short_pathways": "Pathways",
        "clear_chat": "\U0001F5D1\uFE0F  Clear Chat",
        "example_queries": "**Example queries:**",
        "focus_caption": "Focus: {label} | Ask about gene knockouts, brain expression, drug perturbations, receptor binding",
        "chat_placeholder": "Ask about CRISPR knockouts, NDD genes, drug perturbations...",
        "spinner": "Querying {datasets}...",
        "query_log": "### \U0001F4DC Query Log",
        "query_log_caption": "Recent queries from all users",
        "btn_refresh": "\U0001F504 Refresh",
        "btn_clear_log": "\U0001F5D1\uFE0F Clear Log",
        "clear_log_warning": "\u26A0\uFE0F **This action is irreversible.** All query history will be permanently deleted and cannot be recovered.",
        "clear_log_confirm": "Yes, clear log",
        "clear_log_cancel": "Cancel",
        "clear_log_success": "Query log cleared.",
        "clear_log_fail": "Failed to clear log: {err}",
        "no_history": "No queries logged yet. Ask a question to get started!",
        "all_datasets": "All datasets",
        "ds_short_mouse": "Mouse",
        "ds_short_human": "Human CRISPRbrain",
        "ds_short_expression": "Expression",
        "ds_short_pharma": "Pharmacology",
        "ds_short_tahoe": "Tahoe-100M",
        "examples_mouse": [
            "What happens when Psen1 is knocked out?",
            "Which KOs affect Hcrt expression?",
            "Show GABAergic cell types",
            "List knockout targets matching App",
        ],
        "examples_human": [
            "What screens show PSEN1 as a hit?",
            "List human neuron screens",
            "Is SOD2 protective for neuron survival?",
            "What is MAPT expression in brain?",
        ],
        "examples_tahoe": [
            "What drugs are in Tahoe-100M?",
            "Show suvorexant cluster distribution",
            "List orexin antagonists in Tahoe-100M",
            "What MOA classes are represented?",
        ],
        "examples_pharma": [
            "Compare suvorexant vs lemborexant binding",
            "Show OX2R selectivity for DORAs",
            "What genes does lemborexant perturb?",
            "Compare DORA potencies across OX1R and OX2R",
        ],
        "examples_default": [
            "Compare PSEN1 across mouse and human",
            "Summarize all datasets",
            "List ASD risk genes",
            "Is APP important for neuron survival?",
            "What drugs target orexin receptors?",
        ],
    },
    "ja": {
        "subtitle": "\u7a2e\u6a2a\u65ad\u647e\u4e71\u30c7\u30fc\u30bf\u30a8\u30af\u30b9\u30d7\u30ed\u30fc\u30e9\u30fc",
        "sql_ready": "SQL Warehouse: {msg}",
        "sql_error": "SQL: {msg}",
        "dataset_focus": "\U0001F50D \u30c7\u30fc\u30bf\u30bb\u30c3\u30c8\u9078\u629e",
        "btn_all": "\u5168\u9078\u629e",
        "btn_none": "\u5168\u89e3\u9664",
        "ds_mouse": "\U0001F42D \u30de\u30a6\u30b9 (CRISPR Atlas + IMPC)",
        "ds_human": "\U0001F9EC \u30d2\u30c8 CRISPRbrain",
        "ds_expression": "\U0001F4CA \u907a\u4f1d\u5b50\u767a\u73fe (GTEx + \u30b3\u30db\u30fc\u30c8)",
        "ds_pharma": "\U0001F48A \u85ac\u7406\u5b66 (LINCS + ChEMBL)",
        "ds_tahoe": "\U0001F52C Tahoe-100M \u85ac\u5264\u30a2\u30c8\u30e9\u30b9",
        "clear_chat": "\U0001F5D1\uFE0F  \u30c1\u30e3\u30c3\u30c8\u3092\u30af\u30ea\u30a2",
        "example_queries": "**\u8cea\u554f\u4f8b:**",
        "focus_caption": "\u30d5\u30a9\u30fc\u30ab\u30b9: {label} | \u907a\u4f1d\u5b50\u30ce\u30c3\u30af\u30a2\u30a6\u30c8\u3001\u8133\u767a\u73fe\u3001\u85ac\u5264\u647e\u4e71\u3001\u53d7\u5bb9\u4f53\u7d50\u5408\u306b\u3064\u3044\u3066\u8cea\u554f\u3067\u304d\u307e\u3059",
        "chat_placeholder": "CRISPR\u30ce\u30c3\u30af\u30a2\u30a6\u30c8\u3001NDD\u907a\u4f1d\u5b50\u3001\u85ac\u5264\u647e\u4e71\u306b\u3064\u3044\u3066\u8cea\u554f...",
        "spinner": "{datasets} \u3092\u30af\u30a8\u30ea\u4e2d...",
        "query_log": "### \U0001F4DC \u30af\u30a8\u30ea\u30ed\u30b0",
        "query_log_caption": "\u5168\u30e6\u30fc\u30b6\u30fc\u306e\u6700\u8fd1\u306e\u30af\u30a8\u30ea",
        "btn_refresh": "\U0001F504 \u66f4\u65b0",
        "btn_clear_log": "\U0001F5D1\uFE0F \u30ed\u30b0\u524a\u9664",
        "clear_log_warning": "\u26A0\uFE0F **\u3053\u306e\u64cd\u4f5c\u306f\u5143\u306b\u623b\u305b\u307e\u305b\u3093\u3002** \u5168\u3066\u306e\u30af\u30a8\u30ea\u5c65\u6b74\u304c\u5b8c\u5168\u306b\u524a\u9664\u3055\u308c\u3001\u5fa9\u5143\u3067\u304d\u307e\u305b\u3093\u3002",
        "clear_log_confirm": "\u306f\u3044\u3001\u524a\u9664\u3059\u308b",
        "clear_log_cancel": "\u30ad\u30e3\u30f3\u30bb\u30eb",
        "clear_log_success": "\u30af\u30a8\u30ea\u30ed\u30b0\u3092\u524a\u9664\u3057\u307e\u3057\u305f\u3002",
        "clear_log_fail": "\u30ed\u30b0\u524a\u9664\u306b\u5931\u6557\u3057\u307e\u3057\u305f: {err}",
        "no_history": "\u307e\u3060\u30af\u30a8\u30ea\u304c\u3042\u308a\u307e\u305b\u3093\u3002\u8cea\u554f\u3057\u3066\u307f\u307e\u3057\u3087\u3046\uff01",
        "all_datasets": "\u5168\u30c7\u30fc\u30bf\u30bb\u30c3\u30c8",
        "ds_short_mouse": "\u30de\u30a6\u30b9",
        "ds_short_human": "\u30d2\u30c8 CRISPRbrain",
        "ds_short_expression": "\u907a\u4f1d\u5b50\u767a\u73fe",
        "ds_short_pharma": "\u85ac\u7406\u5b66",
        "ds_short_tahoe": "Tahoe-100M",
        "examples_mouse": [
            "Psen1\u3092\u30ce\u30c3\u30af\u30a2\u30a6\u30c8\u3059\u308b\u3068\u4f55\u304c\u8d77\u3053\u308b\uff1f",
            "Hcrt\u767a\u73fe\u306b\u5f71\u97ff\u3059\u308bKO\u306f\uff1f",
            "GABA\u4f5c\u52d5\u6027\u7d30\u80de\u30bf\u30a4\u30d7\u3092\u8868\u793a",
            "App\u306b\u4e00\u81f4\u3059\u308bKO\u30bf\u30fc\u30b2\u30c3\u30c8\u3092\u30ea\u30b9\u30c8",
        ],
        "examples_human": [
            "PSEN1\u304c\u30d2\u30c3\u30c8\u306e\u30b9\u30af\u30ea\u30fc\u30f3\u306f\uff1f",
            "\u30d2\u30c8\u795e\u7d4c\u30b9\u30af\u30ea\u30fc\u30f3\u4e00\u89a7",
            "SOD2\u306f\u795e\u7d4c\u7d30\u80de\u751f\u5b58\u306b\u4fdd\u8b77\u7684\u304b\uff1f",
            "MAPT\u306e\u8133\u5185\u767a\u73fe\u306f\uff1f",
        ],
        "examples_tahoe": [
            "Tahoe-100M\u306b\u542b\u307e\u308c\u308b\u85ac\u5264\u306f\uff1f",
            "\u30b9\u30dc\u30ec\u30ad\u30b5\u30f3\u30c8\u306e\u30af\u30e9\u30b9\u30bf\u30fc\u5206\u5e03\u3092\u8868\u793a",
            "Tahoe-100M\u306e\u30aa\u30ec\u30ad\u30b7\u30f3\u62ee\u6297\u85ac\u3092\u30ea\u30b9\u30c8",
            "\u3069\u306eMOA\u30af\u30e9\u30b9\u304c\u542b\u307e\u308c\u3066\u3044\u308b\u304b\uff1f",
        ],
        "examples_pharma": [
            "\u30b9\u30dc\u30ec\u30ad\u30b5\u30f3\u30c8 vs \u30ec\u30f3\u30dc\u30ec\u30ad\u30b5\u30f3\u30c8\u306e\u7d50\u5408\u6bd4\u8f03",
            "DORA\u306eOX2R\u9078\u629e\u6027\u3092\u8868\u793a",
            "\u30ec\u30f3\u30dc\u30ec\u30ad\u30b5\u30f3\u30c8\u304c\u647e\u4e71\u3059\u308b\u907a\u4f1d\u5b50\u306f\uff1f",
            "OX1R\u3068OX2R\u306b\u304a\u3051\u308bDORA\u306e\u6d3b\u6027\u6bd4\u8f03",
        ],
        "examples_default": [
            "PSEN1\u3092\u30de\u30a6\u30b9\u3068\u30d2\u30c8\u3067\u6bd4\u8f03",
            "\u5168\u30c7\u30fc\u30bf\u30bb\u30c3\u30c8\u306e\u6982\u8981",
            "ASD\u30ea\u30b9\u30af\u907a\u4f1d\u5b50\u3092\u30ea\u30b9\u30c8",
            "APP\u306f\u795e\u7d4c\u7d30\u80de\u751f\u5b58\u306b\u91cd\u8981\u304b\uff1f",
            "\u30aa\u30ec\u30ad\u30b7\u30f3\u53d7\u5bb9\u4f53\u3092\u6a19\u7684\u3068\u3059\u308b\u85ac\u5264\u306f\uff1f",
        ],
    },
}


def _t(key: str, **kwargs) -> str | list:
    """Get translated string for the active language."""
    lang = st.session_state.get("lang", "en")
    val = _I18N.get(lang, _I18N["en"]).get(key, _I18N["en"].get(key, key))
    if isinstance(val, str) and kwargs:
        return val.format(**kwargs)
    return val


# ── SQL Warm-Up ──
@st.cache_resource(show_spinner=False)
def _check_sql_connection():
    try:
        start = time.time()
        _execute_query("SELECT 1 AS ping", timeout=30)
        return True, f"Ready ({time.time() - start:.0f}s)"
    except Exception as e:
        return False, f"{type(e).__name__}: {str(e)[:80]}"

_sql_ok, _sql_msg = _check_sql_connection()

# ── Query Log ──
TABLE_QUERY_LOG = CFG.query_log_fqn


def _get_user_email() -> str:
    """Get current user email from Databricks OAuth context."""
    try:
        from databricks.sdk import WorkspaceClient
        w = WorkspaceClient()
        return w.current_user.me().user_name or "unknown"
    except Exception:
        return "unknown"


def log_query(user_email: str, species_focus: str, query_text: str, response_text: str, tools_used: str = ""):
    """Persist a query and response to the log table."""
    try:
        query_id = str(uuid.uuid4())
        ts = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
        q_esc = query_text.replace("'", "''")
        r_esc = response_text[:4000].replace("'", "''")
        t_esc = tools_used.replace("'", "''")
        sql = f"""
            INSERT INTO {TABLE_QUERY_LOG}
            (query_id, timestamp, user_email, species_focus, query_text, response_text, tools_used)
            VALUES ('{query_id}', '{ts}', '{user_email}', '{species_focus}', '{q_esc}', '{r_esc}', '{t_esc}')
        """
        _execute_query(sql, timeout=15)
    except Exception as e:
        logger.warning(f"Failed to log query: {e}")


@st.cache_data(ttl=30)
def fetch_query_history(limit: int = 50) -> list[dict]:
    """Fetch recent queries from all users."""
    try:
        sql = f"""
            SELECT query_id, timestamp, user_email, species_focus, query_text, response_text
            FROM {TABLE_QUERY_LOG}
            ORDER BY timestamp DESC
            LIMIT {limit}
        """
        return _execute_query(sql, timeout=15)
    except Exception:
        return []


# ── gnomAD + OpenTargets Query Functions ──
_TABLE_GNOMAD = CFG.prefixed_fqn("gnomad")
_TABLE_OPENTARGETS = CFG.prefixed_fqn("opentargets")


def query_gnomad_constraint(gene_symbol: str = "", limit: int = 20) -> str:
    """Query gnomAD gene constraint scores and ClinVar variants."""
    import json as _json
    conditions = ["source_key = 'gnomad'"]
    if gene_symbol:
        conditions.append(f"UPPER(gene_symbol) = UPPER('{gene_symbol}')")
    where = " AND ".join(conditions)
    sql = f"""
        SELECT gene_symbol, title, summary,
               payload:constraint.pLI::DOUBLE AS pLI,
               payload:constraint.oe_lof_upper::DOUBLE AS LOEUF,
               payload:constraint.mis_z::DOUBLE AS mis_z,
               payload:constraint.obs_lof::INT AS obs_lof,
               payload:constraint.exp_lof::DOUBLE AS exp_lof
        FROM {_TABLE_GNOMAD}
        WHERE {where} AND title LIKE '%Gene Constraint%'
        ORDER BY pLI DESC
        LIMIT {limit}
    """
    results = _execute_query(sql, timeout=15)
    return _json.dumps({"source": "gnomAD", "gene": gene_symbol, "count": len(results), "results": results})


def query_gnomad_variants(gene_symbol: str = "", significance: str = "", limit: int = 20) -> str:
    """Query gnomAD ClinVar variants for a gene."""
    import json as _json
    conditions = ["source_key = 'gnomad'", "title NOT LIKE '%Gene Constraint%'"]
    if gene_symbol:
        conditions.append(f"UPPER(gene_symbol) = UPPER('{gene_symbol}')")
    if significance:
        conditions.append(f"LOWER(title) LIKE LOWER('%{significance}%')")
    where = " AND ".join(conditions)
    sql = f"""
        SELECT gene_symbol, title, summary
        FROM {_TABLE_GNOMAD}
        WHERE {where}
        LIMIT {limit}
    """
    results = _execute_query(sql, timeout=15)
    return _json.dumps({"source": "gnomAD ClinVar", "gene": gene_symbol, "count": len(results), "results": results})


def query_opentargets_associations(gene_symbol: str = "", disease: str = "", limit: int = 20) -> str:
    """Query OpenTargets target-disease associations."""
    import json as _json
    conditions = ["source_key = 'opentargets'", "disease IS NOT NULL"]
    if gene_symbol:
        conditions.append(f"UPPER(gene_symbol) = UPPER('{gene_symbol}')")
    if disease:
        conditions.append(f"LOWER(disease) LIKE LOWER('%{disease}%')")
    where = " AND ".join(conditions)
    sql = f"""
        SELECT gene_symbol, disease, title, summary
        FROM {_TABLE_OPENTARGETS}
        WHERE {where}
        ORDER BY title DESC
        LIMIT {limit}
    """
    results = _execute_query(sql, timeout=15)
    return _json.dumps({"source": "OpenTargets", "gene": gene_symbol, "count": len(results), "results": results})


def query_opentargets_tractability(gene_symbol: str = "", limit: int = 20) -> str:
    """Query OpenTargets tractability/druggability for a gene."""
    import json as _json
    conditions = ["source_key = 'opentargets'", "title LIKE '%Tractability%'"]
    if gene_symbol:
        conditions.append(f"UPPER(gene_symbol) = UPPER('{gene_symbol}')")
    where = " AND ".join(conditions)
    sql = f"""
        SELECT gene_symbol, title, summary
        FROM {_TABLE_OPENTARGETS}
        WHERE {where}
        LIMIT {limit}
    """
    results = _execute_query(sql, timeout=15)
    return _json.dumps({"source": "OpenTargets Tractability", "gene": gene_symbol, "count": len(results), "results": results})


# ── New Source Query Functions (CIViC, cBioPortal, GTR, HPA, Monarch, UniProt, Reactome, KEGG) ──
_TABLE_CIVIC = CFG.prefixed_fqn("civic")
_TABLE_CBIOPORTAL = CFG.prefixed_fqn("cbioportal")
_TABLE_NCBI_GTR = CFG.prefixed_fqn("ncbi_gtr")
_TABLE_HPA = CFG.prefixed_fqn("human_protein_atlas")
_TABLE_MONARCH = CFG.prefixed_fqn("monarch")
_TABLE_UNIPROT = CFG.prefixed_fqn("uniprot")
_TABLE_REACTOME = CFG.prefixed_fqn("reactome")
_TABLE_KEGG = CFG.prefixed_fqn("kegg")


def query_civic_evidence(gene_symbol: str = "", disease: str = "", limit: int = 20) -> str:
    """Query CIViC clinical variant evidence and therapy associations."""
    import json as _json
    conditions = ["source_key = 'civic'"]
    if gene_symbol:
        conditions.append(f"UPPER(gene_symbol) = UPPER('{gene_symbol}')")
    if disease:
        conditions.append(f"LOWER(disease) LIKE LOWER('%{disease}%')")
    where = " AND ".join(conditions)
    sql = f"SELECT gene_symbol, disease, drug, title, summary FROM {_TABLE_CIVIC} WHERE {where} LIMIT {limit}"
    results = _execute_query(sql, timeout=15)
    return _json.dumps({"source": "CIViC", "gene": gene_symbol, "count": len(results), "results": results})


def query_cbioportal_mutations(gene_symbol: str = "", cancer_type: str = "", limit: int = 20) -> str:
    """Query cBioPortal cancer mutation data."""
    import json as _json
    conditions = ["source_key = 'cbioportal'"]
    if gene_symbol:
        conditions.append(f"UPPER(gene_symbol) = UPPER('{gene_symbol}')")
    if cancer_type:
        conditions.append(f"LOWER(disease) LIKE LOWER('%{cancer_type}%')")
    where = " AND ".join(conditions)
    sql = f"SELECT gene_symbol, disease, title, summary FROM {_TABLE_CBIOPORTAL} WHERE {where} LIMIT {limit}"
    results = _execute_query(sql, timeout=15)
    return _json.dumps({"source": "cBioPortal", "gene": gene_symbol, "count": len(results), "results": results})


def query_genetic_tests(gene_symbol: str = "", limit: int = 20) -> str:
    """Query NCBI GTR for available genetic tests."""
    import json as _json
    conditions = ["source_key = 'ncbi_gtr'"]
    if gene_symbol:
        conditions.append(f"UPPER(gene_symbol) = UPPER('{gene_symbol}')")
    where = " AND ".join(conditions)
    sql = f"SELECT gene_symbol, title, summary FROM {_TABLE_NCBI_GTR} WHERE {where} LIMIT {limit}"
    results = _execute_query(sql, timeout=15)
    return _json.dumps({"source": "NCBI GTR", "gene": gene_symbol, "count": len(results), "results": results})


def query_protein_expression(gene_symbol: str = "", tissue: str = "", limit: int = 20) -> str:
    """Query Human Protein Atlas tissue expression and localization."""
    import json as _json
    conditions = ["source_key = 'human_protein_atlas'"]
    if gene_symbol:
        conditions.append(f"UPPER(gene_symbol) = UPPER('{gene_symbol}')")
    if tissue:
        conditions.append(f"LOWER(title) LIKE LOWER('%{tissue}%')")
    where = " AND ".join(conditions)
    sql = f"SELECT gene_symbol, title, summary FROM {_TABLE_HPA} WHERE {where} LIMIT {limit}"
    results = _execute_query(sql, timeout=15)
    return _json.dumps({"source": "Human Protein Atlas", "gene": gene_symbol, "count": len(results), "results": results})


def query_monarch_associations(gene_symbol: str = "", disease: str = "", limit: int = 20) -> str:
    """Query Monarch Initiative gene-disease and gene-phenotype associations."""
    import json as _json
    conditions = ["source_key = 'monarch'"]
    if gene_symbol:
        conditions.append(f"UPPER(gene_symbol) = UPPER('{gene_symbol}')")
    if disease:
        conditions.append(f"LOWER(disease) LIKE LOWER('%{disease}%') OR LOWER(title) LIKE LOWER('%{disease}%')")
    where = " AND ".join(conditions)
    sql = f"SELECT gene_symbol, disease, title, summary FROM {_TABLE_MONARCH} WHERE {where} LIMIT {limit}"
    results = _execute_query(sql, timeout=15)
    return _json.dumps({"source": "Monarch Initiative", "gene": gene_symbol, "count": len(results), "results": results})


def query_protein_function(gene_symbol: str = "", limit: int = 10) -> str:
    """Query UniProt protein function, domains, and disease associations."""
    import json as _json
    conditions = ["source_key = 'uniprot'"]
    if gene_symbol:
        conditions.append(f"UPPER(gene_symbol) = UPPER('{gene_symbol}')")
    where = " AND ".join(conditions)
    sql = f"SELECT gene_symbol, disease, title, summary FROM {_TABLE_UNIPROT} WHERE {where} LIMIT {limit}"
    results = _execute_query(sql, timeout=15)
    return _json.dumps({"source": "UniProt", "gene": gene_symbol, "count": len(results), "results": results})


def query_reactome_pathways(gene_symbol: str = "", pathway: str = "", limit: int = 20) -> str:
    """Query Reactome biological pathways for a gene."""
    import json as _json
    conditions = ["source_key = 'reactome'"]
    if gene_symbol:
        conditions.append(f"UPPER(gene_symbol) = UPPER('{gene_symbol}')")
    if pathway:
        conditions.append(f"LOWER(title) LIKE LOWER('%{pathway}%')")
    where = " AND ".join(conditions)
    sql = f"SELECT gene_symbol, title, summary FROM {_TABLE_REACTOME} WHERE {where} LIMIT {limit}"
    results = _execute_query(sql, timeout=15)
    return _json.dumps({"source": "Reactome", "gene": gene_symbol, "count": len(results), "results": results})


def query_kegg_pathways(gene_symbol: str = "", pathway: str = "", limit: int = 20) -> str:
    """Query KEGG metabolic/signaling pathways for a gene."""
    import json as _json
    conditions = ["source_key = 'kegg'"]
    if gene_symbol:
        conditions.append(f"UPPER(gene_symbol) = UPPER('{gene_symbol}')")
    if pathway:
        conditions.append(f"LOWER(title) LIKE LOWER('%{pathway}%')")
    where = " AND ".join(conditions)
    sql = f"SELECT gene_symbol, title, summary FROM {_TABLE_KEGG} WHERE {where} LIMIT {limit}"
    results = _execute_query(sql, timeout=15)
    return _json.dumps({"source": "KEGG", "gene": gene_symbol, "count": len(results), "results": results})


# ── Page Config ──
st.set_page_config(page_title="NeuroPlex", page_icon="favicon.svg", layout="wide")


# ── Tool Registry ──
# ── IMPC MCP Integration ───────────────────────────────────────────────────

_IMPC_BASE = "https://www.ebi.ac.uk/mi/impc/mcp"
_IMPC_HDR  = {"Accept": "application/json, text/event-stream", "Content-Type": "application/json"}


def _call_impc_mcp(endpoint: str, tool_name: str, arguments: dict, timeout: int = 25) -> dict:
    """POST to an IMPC MCP endpoint (SSE JSON-RPC 2.0) and return the parsed result data dict."""
    import requests as _req
    url = f"{_IMPC_BASE}/{endpoint}/"
    payload = {
        "jsonrpc": "2.0", "id": 1,
        "method": "tools/call",
        "params": {"name": tool_name, "arguments": arguments},
    }
    try:
        r = _req.post(url, headers=_IMPC_HDR, json=payload, timeout=timeout, stream=True)
        buf = ""
        for chunk in r.iter_content(chunk_size=None):
            buf += chunk.decode("utf-8")
            if len(buf) > 50000:
                break
        for line in buf.splitlines():
            if line.startswith("data:"):
                rpc = json.loads(line[5:].strip())
                for block in rpc.get("result", {}).get("content", []):
                    if isinstance(block, dict) and block.get("type") == "text":
                        return json.loads(block["text"])
    except Exception as e:
        return {"error": str(e)}
    return {}


def query_impc_phenotypes(gene_symbol: str, phenotype_filter: str = "", limit: int = 15) -> str:
    """Query IMPC for statistically significant in vivo phenotypes in mouse knockouts.
    Returns Mammalian Phenotype (MP) ontology terms, p-values, effect sizes, phenotyping centers."""
    import json as _json
    fq = f"marker_symbol:{gene_symbol}"
    if phenotype_filter:
        fq += f" AND mp_term_name:*{phenotype_filter}*"
    data = _call_impc_mcp("solr", "solr_query", {
        "core": "genotype-phenotype",
        "fq": fq,
        "fl": "marker_symbol,mp_term_name,mp_term_id,p_value,effect_size,phenotyping_center,sex,procedure_name",
        "rows": limit,
        "sort": "p_value asc",
    })
    docs = data.get("response", {}).get("docs", [])
    total = data.get("response", {}).get("numFound", 0)
    return _json.dumps({
        "source": "IMPC genotype-phenotype",
        "gene": gene_symbol,
        "total_phenotypes": total,
        "count": len(docs),
        "results": docs,
    })


def query_impc_orthology(gene_symbol: str, direction: str = "mouse_to_human") -> str:
    """Look up IMPC curated mouse<->human one-to-one ortholog mapping.
    Returns support count, orthology category (GOOD/MODERATE), HGNC and MGI accession IDs."""
    import json as _json
    tool_name = "find_orthologs_by_human_genes" if direction == "human_to_mouse" else "find_orthologs_by_mouse_genes"
    data = _call_impc_mcp("orthology", tool_name, {"genes": gene_symbol})
    return _json.dumps({
        "source": "IMPC orthology",
        "query_gene": gene_symbol,
        "direction": direction,
        "count": len(data.get("data", [])),
        "results": data.get("data", []),
    })


def query_impc_alleles(gene_symbol: str, product_type: str = "mice") -> str:
    """Retrieve available IMPC allele and mouse line products for a gene.
    product_type: mice (live lines, default) or es_cell.
    Returns allele names, types, production centers, and repository ordering links."""
    import json as _json
    tools = {
        "mice":    ("get_mice_products",    {"gene_symbol": gene_symbol}),
        "es_cell": ("get_es_cell_products", {"gene_symbol": gene_symbol}),
    }
    tool_name, args = tools.get(product_type, tools["mice"])
    data = _call_impc_mcp("allele", tool_name, args)
    return _json.dumps({
        "source": "IMPC allele registry",
        "gene": gene_symbol,
        "product_type": product_type,
        "count": len(data.get("data", [])),
        "results": data.get("data", []),
    })


def query_impc_publications(query: str, limit: int = 5) -> str:
    """Semantic search IMPC publications for mouse phenotype studies.
    Finds primary IMPC research articles on specific genes, phenotypes, or disease models."""
    import json as _json
    data = _call_impc_mcp("publication", "by_search_query", {"searchQuery": query})
    results = data.get("data", [])[:limit]
    return _json.dumps({
        "source": "IMPC publications",
        "query": query,
        "count": len(results),
        "results": results,
    })


TOOL_FUNCTIONS = {
    "query_knockout_effects": lambda **kw: query_knockout_effects(**kw),
    "find_knockouts_affecting_gene": lambda **kw: find_knockouts_affecting_gene(**kw),
    "list_ndd_risk_genes": lambda **kw: list_ndd_risk_genes(**kw),
    "list_brain_cell_types": lambda **kw: list_brain_cell_types(**kw),
    "list_knockout_targets": lambda **kw: list_knockout_targets(**kw),
    "query_human_screen_hits": lambda **kw: query_human_screen_hits(**kw),
    "list_human_screens": lambda **kw: list_human_screens(**kw),
    "compare_cross_species": lambda **kw: compare_cross_species(**kw),
    "query_human_expression": lambda **kw: query_human_expression(**kw),
    "query_baseline_expression": lambda **kw: query_baseline_expression(**kw),
    "query_drug_perturbations": lambda **kw: query_drug_perturbations(**kw),
    "query_receptor_pharmacology": lambda **kw: query_receptor_pharmacology(**kw),
    "summarize_atlas": lambda **kw: summarize_atlas(),
    "list_tahoe_drugs": lambda **kw: list_tahoe_drugs(**kw),
    "query_tahoe_drug_clusters": lambda **kw: query_tahoe_drug_clusters(**kw),
    "query_tahoe_moa_distribution": lambda **kw: query_tahoe_moa_distribution(**kw),
    "query_gnomad_constraint": lambda **kw: query_gnomad_constraint(**kw),
    "query_gnomad_variants": lambda **kw: query_gnomad_variants(**kw),
    "query_opentargets_associations": lambda **kw: query_opentargets_associations(**kw),
    "query_opentargets_tractability": lambda **kw: query_opentargets_tractability(**kw),
    "query_civic_evidence": lambda **kw: query_civic_evidence(**kw),
    "query_cbioportal_mutations": lambda **kw: query_cbioportal_mutations(**kw),
    "query_genetic_tests": lambda **kw: query_genetic_tests(**kw),
    "query_protein_expression": lambda **kw: query_protein_expression(**kw),
    "query_monarch_associations": lambda **kw: query_monarch_associations(**kw),
    "query_protein_function": lambda **kw: query_protein_function(**kw),
    "query_reactome_pathways": lambda **kw: query_reactome_pathways(**kw),
    "query_kegg_pathways": lambda **kw: query_kegg_pathways(**kw),
    "query_impc_phenotypes": lambda **kw: query_impc_phenotypes(**kw),
    "query_impc_orthology": lambda **kw: query_impc_orthology(**kw),
    "query_impc_alleles": lambda **kw: query_impc_alleles(**kw),
    "query_impc_publications": lambda **kw: query_impc_publications(**kw),
}

TOOL_SCHEMAS = [
    {"type": "function", "function": {"name": "query_knockout_effects", "description": "[MOUSE] Find genes affected when a gene is knocked out in mouse brain (in vivo). Returns DE genes with log2FC per cell type.", "parameters": {"type": "object", "properties": {"gene_target": {"type": "string", "description": "Mouse gene (Title Case: Psen1, App, Hcrt)"}, "cell_type": {"type": "string", "description": "Cell type filter (Glut, GABA, CTX)"}, "pval_threshold": {"type": "number"}, "limit": {"type": "integer"}}, "required": ["gene_target"]}}},
    {"type": "function", "function": {"name": "find_knockouts_affecting_gene", "description": "[MOUSE] Which knockouts affect a specific gene in mouse brain? Finds upstream regulators.", "parameters": {"type": "object", "properties": {"gene_name": {"type": "string", "description": "Mouse gene (Title Case: Hcrt, Grin1)"}, "cell_type": {"type": "string"}, "pval_threshold": {"type": "number"}, "limit": {"type": "integer"}}, "required": ["gene_name"]}}},
    {"type": "function", "function": {"name": "list_ndd_risk_genes", "description": "NDD/ASD/DD risk genes from TADA analysis.", "parameters": {"type": "object", "properties": {"disorder": {"type": "string", "description": "ASD, DD, or NDD"}, "fdr_threshold": {"type": "number"}, "limit": {"type": "integer"}}, "required": []}}},
    {"type": "function", "function": {"name": "list_brain_cell_types", "description": "[MOUSE] Brain cell types and regions in CRISPR Atlas with cell counts.", "parameters": {"type": "object", "properties": {"region_filter": {"type": "string"}, "limit": {"type": "integer"}}, "required": []}}},
    {"type": "function", "function": {"name": "list_knockout_targets", "description": "[MOUSE] List available CRISPR knockout gene targets.", "parameters": {"type": "object", "properties": {"search": {"type": "string"}, "limit": {"type": "integer"}}, "required": []}}},
    {"type": "function", "function": {"name": "query_human_screen_hits", "description": "[HUMAN] Query CRISPRbrain screen results for a gene in iPSC-derived neurons/microglia/astrocytes. Returns phenotype scores and hit class.", "parameters": {"type": "object", "properties": {"gene": {"type": "string", "description": "Human gene (UPPER CASE: PSEN1, APP, SOD2)"}, "cell_type": {"type": "string", "description": "Cell type (Glutamatergic, Microglia, Astrocyte)"}, "phenotype": {"type": "string", "description": "Phenotype (Survival, Tau, TDP-43)"}, "limit": {"type": "integer"}}, "required": ["gene"]}}},
    {"type": "function", "function": {"name": "list_human_screens", "description": "[HUMAN] List available CRISPRbrain screens with hit counts.", "parameters": {"type": "object", "properties": {"cell_type": {"type": "string"}, "limit": {"type": "integer"}}, "required": []}}},
    {"type": "function", "function": {"name": "compare_cross_species", "description": "[CROSS-SPECIES] Compare a gene across MOUSE (CRISPR Atlas) and HUMAN (CRISPRbrain). Shows conservation of function.", "parameters": {"type": "object", "properties": {"gene": {"type": "string", "description": "Gene symbol (auto-converts case per species)"}, "limit": {"type": "integer"}}, "required": ["gene"]}}},
    {"type": "function", "function": {"name": "query_human_expression", "description": "[HUMAN] Gene expression from BioFINDER/ROSMAP neurology cohorts (brain + blood).", "parameters": {"type": "object", "properties": {"gene_symbol": {"type": "string", "description": "Human gene (HCRT, APP, MAPT)"}, "cohort": {"type": "string"}, "tissue": {"type": "string"}}, "required": ["gene_symbol"]}}},
    {"type": "function", "function": {"name": "query_baseline_expression", "description": "[HUMAN] Baseline gene expression across brain regions from GTEx v8. Shows median TPM per brain subregion. Key for understanding where HCRT/HCRTR1/HCRTR2 are expressed.", "parameters": {"type": "object", "properties": {"gene_symbol": {"type": "string", "description": "Gene symbol (UPPER CASE: HCRT, HCRTR1, HCRTR2, APP)"}, "tissue": {"type": "string", "description": "Brain region filter (Hypothalamus, Cortex, Cerebellum)"}}, "required": ["gene_symbol"]}}},
    {"type": "function", "function": {"name": "query_drug_perturbations", "description": "[PHARMACOLOGY] LINCS L1000 drug perturbation signatures. Shows how orexin antagonists (suvorexant, lemborexant, daridorexant) alter gene expression in neuronal cell lines.", "parameters": {"type": "object", "properties": {"compound": {"type": "string", "description": "Drug name (suvorexant, lemborexant, daridorexant, SB-334867, almorexant)"}, "gene_symbol": {"type": "string", "description": "Gene to check effect on (HCRT, GABRA1, BDNF)"}, "cell_line": {"type": "string", "description": "Cell line (SHSY5Y, NPC, NEU)"}, "neuronal_only": {"type": "boolean", "description": "Restrict to neuronal lines only"}, "limit": {"type": "integer"}}, "required": []}}},
    {"type": "function", "function": {"name": "query_receptor_pharmacology", "description": "[PHARMACOLOGY] ChEMBL binding/activity data for orexin receptor ligands. Returns IC50, Ki, Kd values for compounds vs OX1R/OX2R. Covers suvorexant, lemborexant, daridorexant, seltorexant.", "parameters": {"type": "object", "properties": {"compound": {"type": "string", "description": "Drug name (suvorexant, lemborexant, daridorexant)"}, "target": {"type": "string", "description": "Receptor (OX1R, OX2R, HCRTR1, HCRTR2)"}, "assay_type": {"type": "string", "description": "Activity type (Ki, IC50, Kd)"}, "limit": {"type": "integer"}}, "required": []}}},
    {"type": "function", "function": {"name": "summarize_atlas", "description": "Summary of ALL datasets: Mouse CRISPR Atlas + Human CRISPRbrain + GTEx + LINCS + ChEMBL.", "parameters": {"type": "object", "properties": {}, "required": []}}},
    {"type": "function", "function": {"name": "list_tahoe_drugs", "description": "[TAHOE-100M] List compounds in the Tahoe-100M drug perturbation atlas with MOA and cell counts across cell lines.", "parameters": {"type": "object", "properties": {"search": {"type": "string", "description": "Drug name substring"}, "cell_line": {"type": "string", "description": "Cell line filter"}, "limit": {"type": "integer"}}, "required": []}}},
    {"type": "function", "function": {"name": "query_tahoe_drug_clusters", "description": "[TAHOE-100M] Query transcriptomic cluster distribution for a drug. Shows how drug-treated cells distribute across clusters per cell line.", "parameters": {"type": "object", "properties": {"drug": {"type": "string", "description": "Drug name to query"}, "cell_line": {"type": "string", "description": "Optional cell line filter"}, "limit": {"type": "integer"}}, "required": []}}},
    {"type": "function", "function": {"name": "query_tahoe_moa_distribution", "description": "[TAHOE-100M] Summarize Tahoe-100M compounds by mechanism of action (MOA). Shows drug and cell counts per MOA class.", "parameters": {"type": "object", "properties": {"moa_filter": {"type": "string", "description": "MOA class substring filter"}, "limit": {"type": "integer"}}, "required": []}}},
    {"type": "function", "function": {"name": "query_gnomad_constraint", "description": "[GENETICS] Query gnomAD gene constraint scores (pLI, LOEUF, missense Z). Shows how intolerant a gene is to loss-of-function variants. High pLI (>0.9) = haploinsufficient/essential.", "parameters": {"type": "object", "properties": {"gene_symbol": {"type": "string", "description": "Gene symbol (HCRT, PSEN1, APP)"}, "limit": {"type": "integer"}}, "required": ["gene_symbol"]}}},
    {"type": "function", "function": {"name": "query_gnomad_variants", "description": "[GENETICS] Query gnomAD ClinVar variants for a gene. Returns pathogenic/benign classifications with population allele frequencies.", "parameters": {"type": "object", "properties": {"gene_symbol": {"type": "string", "description": "Gene symbol (PSEN1, APP, MAPT)"}, "significance": {"type": "string", "description": "Clinical significance filter (Pathogenic, Benign, Uncertain)"}, "limit": {"type": "integer"}}, "required": ["gene_symbol"]}}},
    {"type": "function", "function": {"name": "query_opentargets_associations", "description": "[DRUGGABILITY] Query OpenTargets target-disease association scores. Shows evidence-based links between genes and diseases with therapeutic area context.", "parameters": {"type": "object", "properties": {"gene_symbol": {"type": "string", "description": "Gene symbol (HCRT, PSEN1, SOD1)"}, "disease": {"type": "string", "description": "Disease name filter (Alzheimer, narcolepsy, ALS)"}, "limit": {"type": "integer"}}, "required": []}}},
    {"type": "function", "function": {"name": "query_opentargets_tractability", "description": "[DRUGGABILITY] Query OpenTargets druggability/tractability for a gene target. Shows which modalities (small molecule, antibody, PROTAC) are feasible.", "parameters": {"type": "object", "properties": {"gene_symbol": {"type": "string", "description": "Gene symbol (HCRTR1, HCRTR2, PSEN1)"}, "limit": {"type": "integer"}}, "required": ["gene_symbol"]}}},
    {"type": "function", "function": {"name": "query_civic_evidence", "description": "[GENETICS] Query CIViC clinical variant evidence. Shows curated therapy associations, evidence levels (A-E), and clinical significance for gene variants.", "parameters": {"type": "object", "properties": {"gene_symbol": {"type": "string", "description": "Gene symbol (PSEN1, SOD1, SCN2A)"}, "disease": {"type": "string", "description": "Disease filter (Alzheimer, ALS, cancer)"}, "limit": {"type": "integer"}}, "required": ["gene_symbol"]}}},
    {"type": "function", "function": {"name": "query_cbioportal_mutations", "description": "[GENETICS] Query cBioPortal cancer mutations for a gene. Shows mutation types, protein changes, and study context.", "parameters": {"type": "object", "properties": {"gene_symbol": {"type": "string", "description": "Gene symbol (PSEN1, APP, TP53)"}, "cancer_type": {"type": "string", "description": "Cancer type filter (glioblastoma, breast)"}, "limit": {"type": "integer"}}, "required": ["gene_symbol"]}}},
    {"type": "function", "function": {"name": "query_genetic_tests", "description": "[GENETICS] Query NCBI GTR for available clinical genetic tests for a gene.", "parameters": {"type": "object", "properties": {"gene_symbol": {"type": "string", "description": "Gene symbol (PSEN1, MAPT, C9orf72)"}, "limit": {"type": "integer"}}, "required": ["gene_symbol"]}}},
    {"type": "function", "function": {"name": "query_protein_expression", "description": "[EXPRESSION] Query Human Protein Atlas for tissue-level protein expression and subcellular localization.", "parameters": {"type": "object", "properties": {"gene_symbol": {"type": "string", "description": "Gene symbol (HCRT, PSEN1, APP)"}, "tissue": {"type": "string", "description": "Tissue filter (brain, liver, kidney)"}, "limit": {"type": "integer"}}, "required": ["gene_symbol"]}}},
    {"type": "function", "function": {"name": "query_monarch_associations", "description": "[EXPRESSION] Query Monarch Initiative gene-disease associations and phenotype matching from model organisms.", "parameters": {"type": "object", "properties": {"gene_symbol": {"type": "string", "description": "Gene symbol (HCRT, PSEN1, SOD1)"}, "disease": {"type": "string", "description": "Disease or phenotype filter"}, "limit": {"type": "integer"}}, "required": ["gene_symbol"]}}},
    {"type": "function", "function": {"name": "query_protein_function", "description": "[PATHWAYS] Query UniProt protein function, structure, domains, subcellular location, and disease associations.", "parameters": {"type": "object", "properties": {"gene_symbol": {"type": "string", "description": "Gene symbol (HCRT, PSEN1, APP)"}, "limit": {"type": "integer"}}, "required": ["gene_symbol"]}}},
    {"type": "function", "function": {"name": "query_reactome_pathways", "description": "[PATHWAYS] Query Reactome biological pathways a gene participates in. Shows pathway hierarchy, compartments, and reactions.", "parameters": {"type": "object", "properties": {"gene_symbol": {"type": "string", "description": "Gene symbol (HCRT, PSEN1, MAPT)"}, "pathway": {"type": "string", "description": "Pathway name filter (signaling, apoptosis, metabolism)"}, "limit": {"type": "integer"}}, "required": ["gene_symbol"]}}},
    {"type": "function", "function": {"name": "query_kegg_pathways", "description": "[PATHWAYS] Query KEGG metabolic and signaling pathway maps for a gene. Shows pathway class, compounds, and interactions.", "parameters": {"type": "object", "properties": {"gene_symbol": {"type": "string", "description": "Gene symbol (HCRT, PSEN1, SOD1)"}, "pathway": {"type": "string", "description": "Pathway name filter (Alzheimer, orexin, apoptosis)"}, "limit": {"type": "integer"}}, "required": ["gene_symbol"]}}},
    {"type": "function", "function": {"name": "query_impc_phenotypes", "description": "[MOUSE/IMPC] Query the International Mouse Phenotyping Consortium (IMPC) for in vivo phenotypes of mouse gene knockouts. Returns Mammalian Phenotype (MP) ontology terms with p-values, effect sizes, sex, and phenotyping centers from systematic whole-organism screens (locomotor, cardiac, metabolic, neurological).", "parameters": {"type": "object", "properties": {"gene_symbol": {"type": "string", "description": "Mouse gene (Title Case: Hcrt, Psen1, App)"}, "phenotype_filter": {"type": "string", "description": "MP term substring filter (e.g. locomotor, cardiovascular, body weight, brain)"}, "limit": {"type": "integer"}}, "required": ["gene_symbol"]}}},
    {"type": "function", "function": {"name": "query_impc_orthology", "description": "[MOUSE/IMPC] Look up IMPC curated mouse<->human one-to-one ortholog mapping. Returns support count, orthology category (GOOD/MODERATE/WEAK), HGNC and MGI accession IDs. Complements compare_cross_species for validating gene conservation.", "parameters": {"type": "object", "properties": {"gene_symbol": {"type": "string", "description": "Gene symbol: mouse (Title Case: Hcrt, Psen1) or human (UPPER: HCRT, PSEN1)"}, "direction": {"type": "string", "description": "mouse_to_human (default) or human_to_mouse"}}, "required": ["gene_symbol"]}}},
    {"type": "function", "function": {"name": "query_impc_alleles", "description": "[MOUSE/IMPC] Retrieve available IMPC mouse allele and line products for a gene, including live mouse colonies (C57BL/6 background), deletion/reporter alleles, production centers, and ordering links from MMRRC, JAX, and other repositories.", "parameters": {"type": "object", "properties": {"gene_symbol": {"type": "string", "description": "Mouse gene (Title Case: Hcrt, Psen1, App)"}, "product_type": {"type": "string", "description": "mice (live mouse lines, default) or es_cell"}}, "required": ["gene_symbol"]}}},
    {"type": "function", "function": {"name": "query_impc_publications", "description": "[MOUSE/IMPC] Semantic search the IMPC publications database for mouse phenotype studies from IMPC member institutes. Finds primary research articles on specific genes, phenotypes, or disease models.", "parameters": {"type": "object", "properties": {"query": {"type": "string", "description": "Search query (e.g. 'Hcrt sleep narcolepsy', 'Psen1 Alzheimer mouse model')"}, "limit": {"type": "integer"}}, "required": ["query"]}}},
]

SYSTEM_PROMPT = """You are a neuroscience research assistant with access to CRISPR perturbation data,
gene expression, pharmacology, genetics, and druggability datasets spanning EIGHT data sources:

## MOUSE DATA (PerturbAI CRISPR Atlas + IMPC)

### PerturbAI Whole Brain CRISPR Atlas
- 7.7 million single cells from mouse brain (in vivo)
- ~2,046 gene knockouts mapped across hundreds of neuronal cell types
- 745 million pre-computed differential expression results
- Gene symbols: Title Case (Psen1, App, Hcrt, Grin1)
- Tools: query_knockout_effects, find_knockouts_affecting_gene, list_brain_cell_types, list_knockout_targets

### IMPC (International Mouse Phenotyping Consortium)
- Systematic whole-organism knockout phenotyping across hundreds of genes
- Mammalian Phenotype (MP) ontology: locomotor, cardiovascular, metabolic, neurological, body composition
- Orthology: curated mouse<->human one-to-one gene mappings with support scores (GOOD/MODERATE/WEAK)
- Allele registry: live mouse lines (C57BL/6), CRISPR alleles, ES cell clones with repository ordering links
- Publications: semantic search of primary IMPC phenotyping studies
- Tools: query_impc_phenotypes, query_impc_orthology, query_impc_alleles, query_impc_publications

## HUMAN DATA (CRISPRbrain.org)
- iPSC-derived glutamatergic neurons, microglia, astrocytes
- CRISPRi/CRISPRa genome-wide screens
- Phenotypes: survival, tau aggregation, TDP-43, lipid peroxidation, phagocytosis
- Disease models: MAPT-V337M (FTD), PSAP KO
- Gene symbols: UPPER CASE (PSEN1, APP, MAPT, SOD2)
- Tools: query_human_screen_hits, list_human_screens

## HUMAN EXPRESSION (BioFINDER, ROSMAP)
- Full transcriptome RNA-seq from neurology cohorts
- Brain tissue (prefrontal cortex) + blood
- Tool: query_human_expression

## BASELINE BRAIN EXPRESSION (GTEx v8)
- Median TPM expression across 13 brain subregions
- Key for localizing orexin system: HCRT (hypothalamus-specific), HCRTR1/HCRTR2 (widespread)
- Covers AD/NDD risk genes, sleep/circadian genes
- Tool: query_baseline_expression

## DRUG PERTURBATION SIGNATURES (LINCS L1000)
- Transcriptomic effects of orexin receptor antagonists in neuronal cell lines
- Compounds: suvorexant, lemborexant (Eisai/Dayvigo), daridorexant, almorexant, SB-334867, TCS-OX2-29
- Cell lines: SH-SY5Y, NPC, NEU (neuronal) + reference lines
- Z-scores across orexin pathway genes, GABA, glutamate, BDNF signaling
- Tool: query_drug_perturbations

## RECEPTOR PHARMACOLOGY (ChEMBL)
- Binding affinity data (Ki, IC50, Kd, pChEMBL) for orexin receptor ligands
- Targets: OX1R (HCRTR1) and OX2R (HCRTR2)
- Clinical DORAs: suvorexant (Belsomra), lemborexant (Dayvigo/Eisai), daridorexant (Quviviq)
- Selective antagonists: seltorexant (OX2R), SB-334867 (OX1R)
- Tool: query_receptor_pharmacology

## TAHOE-100M DRUG PERTURBATION ATLAS
- ~100 million single cells treated with diverse drug compounds
- Multiple cancer/neuronal cell lines
- MOA (mechanism of action) annotations per compound
- Cluster assignments capturing transcriptomic response states
- SMILES structures + PubChem CIDs available
- Tools: list_tahoe_drugs, query_tahoe_drug_clusters, query_tahoe_moa_distribution

## CROSS-SPECIES
- compare_cross_species: Queries BOTH mouse and human CRISPR data for the same gene
- NDD risk genes: list_ndd_risk_genes (TADA analysis for ASD/DD/NDD/SCZ)

IMPORTANT:
- Always specify species context when presenting results
- Mouse genes: Title Case (Psen1). Human genes: UPPER CASE (PSEN1).
- When user asks about a gene without specifying species, use compare_cross_species
- For orexin/sleep questions, combine baseline expression + pharmacology + CRISPR data
- Lemborexant is Eisai's compound (Dayvigo) — highlight comparative pharmacology when relevant
- The user can select one or more dataset focuses (Mouse/Human/Expression/Pharmacology/Tahoe-100M) from checkboxes in the sidebar
- For Tahoe-100M questions, use list_tahoe_drugs, query_tahoe_drug_clusters, query_tahoe_moa_distribution
- Explain biological significance, cross-species implications, and therapeutic relevance

VISUALIZATION:
- Charts are AUTOMATICALLY rendered below your response when tools return numeric data
- Supported chart types: brain expression bar plots, receptor binding grouped bars, perturbation z-score bars, knockout log2FC bars, cross-species comparisons
- When users ask for visualizations, plots, or comparisons, call the relevant tools — charts appear automatically
- You do NOT need to describe the chart in detail — just explain the biological interpretation
- For best visualizations, query specific genes/compounds rather than broad lists
"""


def get_client():
    from databricks.sdk import WorkspaceClient
    from openai import OpenAI
    w = WorkspaceClient()
    headers = w.config.authenticate()
    token = headers.get("Authorization", "").replace("Bearer ", "")
    host = w.config.host.rstrip("/")
    return OpenAI(api_key=token, base_url=f"{host}/serving-endpoints")


# ── Tool → Dataset Mapping (for hard filtering) ──
_TOOL_DATASET_MAP = {
    "query_knockout_effects": "mouse",
    "find_knockouts_affecting_gene": "mouse",
    "list_brain_cell_types": "mouse",
    "list_knockout_targets": "mouse",
    "query_impc_phenotypes": "mouse",
    "query_impc_orthology": "mouse",
    "query_impc_alleles": "mouse",
    "query_impc_publications": "mouse",
    "query_human_screen_hits": "human",
    "list_human_screens": "human",
    "query_human_expression": "expression",
    "query_baseline_expression": "expression",
    "query_drug_perturbations": "pharma",
    "query_receptor_pharmacology": "pharma",
    "list_tahoe_drugs": "tahoe",
    "query_tahoe_drug_clusters": "tahoe",
    "query_tahoe_moa_distribution": "tahoe",
    "query_gnomad_constraint": "genetics",
    "query_gnomad_variants": "genetics",
    "query_civic_evidence": "genetics",
    "query_cbioportal_mutations": "genetics",
    "query_genetic_tests": "genetics",
    "query_opentargets_associations": "druggability",
    "query_opentargets_tractability": "druggability",
    "query_protein_expression": "protein",
    "query_monarch_associations": "protein",
    "query_protein_function": "pathways",
    "query_reactome_pathways": "pathways",
    "query_kegg_pathways": "pathways",
    # Cross-dataset tools — always available
    "compare_cross_species": None,
    "list_ndd_risk_genes": None,
    "summarize_atlas": None,
}


def run_agent(messages: list[dict], active_datasets: list[str] | None = None, max_iterations: int = 10) -> str:
    client = get_client()
    # Reset tool results for visualization
    st.session_state.tool_results = []

    _ALL_DS = {"mouse", "human", "expression", "pharma", "tahoe", "genetics", "druggability", "protein", "pathways"}
    if not active_datasets:
        active_datasets = list(_ALL_DS)

    active_set = set(active_datasets)

    # Hard-filter tools: only present tools whose dataset is selected (or dataset-agnostic)
    filtered_tools = [
        schema for schema in TOOL_SCHEMAS
        if _TOOL_DATASET_MAP.get(schema["function"]["name"]) is None
        or _TOOL_DATASET_MAP.get(schema["function"]["name"]) in active_set
    ]

    # Build focused system prompt note
    species_note = ""
    if active_set != _ALL_DS:
        _ds_labels = {
            "mouse": "Mouse data: PerturbAI CRISPR Atlas (7.7M cells, 2,046 KOs, 745M DE results) + IMPC whole-organism phenotyping (MP terms, orthology, allele registry, publications)",
            "human": "Human CRISPRbrain — iPSC neurons, microglia, astrocytes; genome-wide CRISPRi/a; phenotypes: survival, tau, TDP-43, phagocytosis",
            "expression": "Brain expression (GTEx baseline + cohort data)",
            "pharma": "Pharmacology (LINCS L1000 + ChEMBL receptor binding). Compare DORAs; highlight Eisai's lemborexant.",
            "tahoe": "Tahoe-100M drug perturbation atlas",
            "genetics": "Genetics (gnomAD gene constraint, pLI/LOEUF scores, ClinVar variant pathogenicity, population allele frequencies)",
            "druggability": "Druggability (OpenTargets target-disease associations, tractability/druggability modalities)",
            "protein": "Protein & Phenotype (Human Protein Atlas tissue expression, Monarch Initiative gene-disease/phenotype associations)",
            "pathways": "Pathways & Function (UniProt protein function/domains, Reactome biological pathways, KEGG metabolic/signaling maps)",
        }
        _ds_tool_hints = {
            "mouse": "Mouse: PerturbAI CRISPR Atlas (in vivo KO effects, cell types, upstream regulators) + IMPC phenotyping (MP terms, orthology, allele lines, publications)",
            "human": "Human CRISPRbrain (iPSC screen hits, phenotype scores)",
            "expression": "GTEx baseline brain expression and BioFINDER/ROSMAP cohort expression",
            "pharma": "Pharmacology — LINCS L1000 drug signatures and ChEMBL receptor binding (Ki/IC50)",
            "tahoe": "Tahoe-100M drug perturbation atlas (cluster distributions, MOA)",
            "genetics": "Genetics — gnomAD gene constraint (pLI, LOEUF, missense Z) and ClinVar variant classifications",
            "druggability": "Druggability — OpenTargets target-disease scores and tractability (small molecule, antibody, PROTAC)",
        }
        active_notes = [_ds_labels[k] for k in active_datasets if k in _ds_labels]
        disabled_keys = _ALL_DS - active_set
        disabled_notes = [_ds_tool_hints[k] for k in sorted(disabled_keys) if k in _ds_tool_hints]
        if active_notes:
            species_note = "\n\nACTIVE DATASETS (only these are available):\n" + "\n".join(f"- {n}" for n in active_notes)
        if disabled_notes:
            species_note += "\n\nDISABLED DATASETS (you cannot query these):\n" + "\n".join(f"- {n}" for n in disabled_notes)
            species_note += (
                "\n\nIMPORTANT: If the user's question primarily requires a DISABLED dataset, "
                "you MUST begin your response with a clear warning: \"⚠️ The dataset needed to answer this "
                "question ([dataset name]) is currently deselected in Dataset Focus. "
                "Please enable it in the sidebar to query this data.\" "
                "Then briefly explain what the disabled dataset contains and how it would answer their question. "
                "Do NOT fabricate an answer from memory — only present data from tools you can actually call."
            )

    # Language instruction
    lang = st.session_state.get("lang", "en")
    lang_instruction = ""
    if lang == "ja":
        lang_instruction = (
            "\n\nLANGUAGE: Respond entirely in Japanese (日本語で回答してください). "
            "Use standard Japanese scientific terminology. "
            "Gene symbols remain in their original form (HCRT, PSEN1, App). "
            "Drug names may use katakana where standard (スボレキサント, レンボレキサント, ダリドレキサント). "
            "Explain biological significance in Japanese. "
            "When warning about disabled datasets, also write the warning in Japanese."
        )

    full_messages = [{"role": "system", "content": SYSTEM_PROMPT + species_note + lang_instruction}] + messages

    _model_name = os.environ.get("DATABRICKS_SERVING_ENDPOINT", CFG.serving_endpoint)

    # Models that do NOT support the temperature parameter
    _NO_TEMPERATURE_MODELS = {"databricks-claude-sonnet-5", "databricks-claude-opus-4"}
    _supports_temperature = not any(m in _model_name for m in _NO_TEMPERATURE_MODELS)

    def _norm_content(raw) -> str:
        """Normalize message content to str.

        Claude models served via Databricks may return content as a list of
        typed blocks (e.g. [TextBlock(type='text', text='...')] or
        [{'type': 'text', 'text': '...'}]) rather than a plain string.
        Flatten to a single string so callers never receive a list.
        """
        if raw is None:
            return ""
        if isinstance(raw, str):
            return raw
        if isinstance(raw, list):
            return "\n".join(
                (b.text if hasattr(b, "text") else b.get("text", ""))
                for b in raw
                if hasattr(b, "text") or (isinstance(b, dict) and b.get("type") == "text")
            )
        return str(raw)  # fallback for any other unexpected type

    for _ in range(max_iterations):
        create_kwargs = dict(
            model=_model_name,
            messages=full_messages,
            tools=filtered_tools,
            tool_choice="auto",
        )
        if _supports_temperature:
            create_kwargs["temperature"] = 0
        response = client.chat.completions.create(**create_kwargs)
        msg = response.choices[0].message
        if not msg.tool_calls:
            return _norm_content(msg.content)

        full_messages.append({
            "role": "assistant",
            "content": _norm_content(msg.content),
            "tool_calls": [{"id": tc.id, "type": "function", "function": {"name": tc.function.name, "arguments": tc.function.arguments}} for tc in msg.tool_calls],
        })

        for tc in msg.tool_calls:
            fn_name = tc.function.name
            fn_args = json.loads(tc.function.arguments)
            fn = TOOL_FUNCTIONS.get(fn_name)
            if fn:
                try:
                    result = fn(**fn_args)
                except Exception as e:
                    result = json.dumps({"error": str(e)})
            else:
                result = json.dumps({"error": f"Unknown tool: {fn_name}"})
            full_messages.append({"role": "tool", "tool_call_id": tc.id, "content": result})
            # Store for visualization
            st.session_state.tool_results.append({"tool": fn_name, "args": fn_args, "result": result})

    return "Max iterations reached."


# ── Visualization Engine ──
import plotly.graph_objects as go
import plotly.express as px


def render_charts(tool_results: list[dict]):
    """Render Plotly charts based on tool call results."""
    import pandas as pd
    import numpy as np

    # Collect all expression results for potential multi-gene heatmap
    expression_results = []

    for tr in tool_results:
        tool = tr["tool"]
        try:
            data = json.loads(tr["result"])
        except (json.JSONDecodeError, TypeError):
            continue

        results = data.get("results", [])
        if not results:
            continue

        # ══════════════════════════════════════════════════════════════════
        # 1. GTEx BASELINE EXPRESSION - Bar + Radar
        # ══════════════════════════════════════════════════════════════════
        if tool == "query_baseline_expression" and len(results) > 1:
            gene = data.get("gene", "Gene")
            tissues = [r.get("tissue", "").replace("Brain - ", "") for r in results]
            tpms = [float(r.get("median_tpm", 0)) for r in results]

            # Store for multi-gene heatmap
            expression_results.append({"gene": gene, "tissues": tissues, "tpms": tpms})

            # Horizontal bar chart
            fig = go.Figure(go.Bar(
                x=tpms, y=tissues, orientation="h",
                marker_color="#a855f7",
                text=[f"{t:.2f}" for t in tpms],
                textposition="outside",
            ))
            fig.update_layout(
                title=f"{gene} Expression Across Brain Regions (GTEx v8)",
                xaxis_title="Median TPM",
                yaxis=dict(autorange="reversed"),
                height=max(300, len(tissues) * 28),
                margin=dict(l=200, r=50, t=50, b=40),
                template="plotly_dark",
            )
            st.plotly_chart(fig, use_container_width=True)

            # Radar chart (brain region profile) - only if >4 regions
            if len(tissues) > 4:
                fig_radar = go.Figure(go.Scatterpolar(
                    r=tpms + [tpms[0]],  # close the polygon
                    theta=tissues + [tissues[0]],
                    fill="toself",
                    fillcolor="rgba(168, 85, 247, 0.3)",
                    line=dict(color="#a855f7", width=2),
                    name=gene,
                ))
                fig_radar.update_layout(
                    title=f"{gene} Brain Region Profile",
                    polar=dict(radialaxis=dict(visible=True, range=[0, max(tpms) * 1.1])),
                    height=400, template="plotly_dark",
                    showlegend=False,
                )
                st.plotly_chart(fig_radar, use_container_width=True)

        # ══════════════════════════════════════════════════════════════════
        # 2. RECEPTOR PHARMACOLOGY - Grouped Bar + Selectivity Scatter + Lollipop
        # ══════════════════════════════════════════════════════════════════
        elif tool == "query_receptor_pharmacology" and len(results) > 1:
            valid = [r for r in results if r.get("pchembl_value")]
            if len(valid) >= 2:
                compounds = [r.get("molecule_name", "") for r in valid]
                targets = [r.get("target_alias", "") for r in valid]
                values = [float(r.get("standard_value", 0)) for r in valid]
                pchembl = [float(r.get("pchembl_value", 0)) for r in valid]
                units = valid[0].get("standard_units", "nM") if valid else "nM"
                assay = valid[0].get("standard_type", "Ki") if valid else "Ki"
                df_viz = pd.DataFrame({"Compound": compounds, "Receptor": targets, f"{assay} ({units})": values, "pChEMBL": pchembl})

                # Grouped bar chart
                fig = px.bar(
                    df_viz, x="Compound", y=f"{assay} ({units})", color="Receptor",
                    barmode="group",
                    color_discrete_map={"OX1R": "#06b6d4", "OX2R": "#fb923c", "HCRTR1": "#06b6d4", "HCRTR2": "#fb923c"},
                )
                fig.update_layout(
                    title=f"Orexin Receptor Binding ({assay})",
                    yaxis_title=f"{assay} ({units})", yaxis_type="log",
                    height=400, template="plotly_dark",
                )
                st.plotly_chart(fig, use_container_width=True)

                # Selectivity scatter (OX1R Ki vs OX2R Ki)
                ox1r = df_viz[df_viz["Receptor"].isin(["OX1R", "HCRTR1"])].groupby("Compound").first().reset_index()
                ox2r = df_viz[df_viz["Receptor"].isin(["OX2R", "HCRTR2"])].groupby("Compound").first().reset_index()
                if len(ox1r) >= 2 and len(ox2r) >= 2:
                    merged = ox1r[["Compound", f"{assay} ({units})"]].merge(
                        ox2r[["Compound", f"{assay} ({units})"]],
                        on="Compound", suffixes=("_OX1R", "_OX2R")
                    )
                    if len(merged) >= 2:
                        fig_scatter = go.Figure()
                        fig_scatter.add_trace(go.Scatter(
                            x=merged[f"{assay} ({units})_OX1R"],
                            y=merged[f"{assay} ({units})_OX2R"],
                            mode="markers+text",
                            text=merged["Compound"],
                            textposition="top center",
                            marker=dict(size=14, color="#a855f7", line=dict(width=2, color="white")),
                        ))
                        # Diagonal line (equipotent)
                        max_val = max(merged[f"{assay} ({units})_OX1R"].max(), merged[f"{assay} ({units})_OX2R"].max()) * 1.5
                        fig_scatter.add_trace(go.Scatter(
                            x=[0.01, max_val], y=[0.01, max_val],
                            mode="lines", line=dict(dash="dash", color="gray", width=1),
                            showlegend=False,
                        ))
                        fig_scatter.update_layout(
                            title=f"OX1R vs OX2R Selectivity ({assay})",
                            xaxis_title=f"OX1R {assay} ({units})", yaxis_title=f"OX2R {assay} ({units})",
                            xaxis_type="log", yaxis_type="log",
                            height=450, template="plotly_dark",
                            annotations=[dict(x=0.9, y=0.1, xref="paper", yref="paper",
                                text="OX2R-selective (below line)", showarrow=False, font=dict(color="#fb923c", size=10)),
                                dict(x=0.1, y=0.9, xref="paper", yref="paper",
                                text="OX1R-selective (above line)", showarrow=False, font=dict(color="#06b6d4", size=10))],
                        )
                        st.plotly_chart(fig_scatter, use_container_width=True)

                # pChEMBL lollipop chart
                if len(df_viz) >= 3:
                    df_lollipop = df_viz.sort_values("pChEMBL", ascending=True)
                    labels = [f"{r['Compound']} ({r['Receptor']})" for _, r in df_lollipop.iterrows()]
                    fig_lollipop = go.Figure()
                    fig_lollipop.add_trace(go.Scatter(
                        x=df_lollipop["pChEMBL"], y=labels,
                        mode="markers",
                        marker=dict(size=12, color=["#06b6d4" if t in ["OX1R", "HCRTR1"] else "#fb923c" for t in df_lollipop["Receptor"]]),
                    ))
                    for i, row in df_lollipop.iterrows():
                        fig_lollipop.add_shape(type="line", x0=0, x1=row["pChEMBL"],
                            y0=labels[list(df_lollipop.index).index(i)], y1=labels[list(df_lollipop.index).index(i)],
                            line=dict(color="gray", width=1))
                    fig_lollipop.update_layout(
                        title="Compound Potency (pChEMBL, higher = more potent)",
                        xaxis_title="pChEMBL value",
                        height=max(300, len(labels) * 30),
                        margin=dict(l=200, r=50, t=50, b=40),
                        template="plotly_dark", showlegend=False,
                    )
                    st.plotly_chart(fig_lollipop, use_container_width=True)

        # ══════════════════════════════════════════════════════════════════
        # 3. LINCS DRUG PERTURBATIONS - Bar + Heatmap
        # ══════════════════════════════════════════════════════════════════
        elif tool == "query_drug_perturbations" and len(results) > 3:
            df_viz = pd.DataFrame(results)
            if "zscore" in df_viz.columns and "gene_symbol" in df_viz.columns:
                df_viz["zscore"] = df_viz["zscore"].astype(float)
                df_top = df_viz.reindex(df_viz["zscore"].abs().sort_values(ascending=False).index[:20])
                df_top = df_top.sort_values("zscore")
                compound = data.get("compound_filter", "Compound")

                # Diverging bar chart
                fig = go.Figure(go.Bar(
                    x=df_top["zscore"], y=df_top["gene_symbol"], orientation="h",
                    marker_color=["#ef4444" if z < 0 else "#22c55e" for z in df_top["zscore"]],
                    text=[f"{z:.2f}" for z in df_top["zscore"]], textposition="outside",
                ))
                fig.update_layout(
                    title=f"Top Gene Perturbations \u2014 {compound}",
                    xaxis_title="Z-score",
                    height=max(350, len(df_top) * 25),
                    margin=dict(l=100, r=50, t=50, b=40),
                    template="plotly_dark",
                )
                st.plotly_chart(fig, use_container_width=True)

                # Heatmap (compound x gene) if multiple compounds present
                if "compound" in df_viz.columns and df_viz["compound"].nunique() > 1:
                    top_genes = df_viz.groupby("gene_symbol")["zscore"].apply(lambda x: x.abs().max()).nlargest(15).index.tolist()
                    df_heat = df_viz[df_viz["gene_symbol"].isin(top_genes)].pivot_table(
                        index="gene_symbol", columns="compound", values="zscore", aggfunc="mean"
                    ).fillna(0)
                    fig_heat = go.Figure(go.Heatmap(
                        z=df_heat.values, x=df_heat.columns.tolist(), y=df_heat.index.tolist(),
                        colorscale="RdBu_r", zmid=0,
                        text=np.round(df_heat.values, 2), texttemplate="%{text}",
                    ))
                    fig_heat.update_layout(
                        title="Drug Perturbation Heatmap (Compound x Gene)",
                        height=max(400, len(top_genes) * 30),
                        margin=dict(l=100, r=50, t=50, b=80),
                        template="plotly_dark",
                    )
                    st.plotly_chart(fig_heat, use_container_width=True)

        # ══════════════════════════════════════════════════════════════════
        # 4. CRISPR KNOCKOUT EFFECTS - Bar + Volcano Plot + Dot Plot
        # ══════════════════════════════════════════════════════════════════
        elif tool == "query_knockout_effects" and len(results) > 3:
            df_viz = pd.DataFrame(results)
            if "log2fc" in df_viz.columns:
                df_viz["log2fc"] = df_viz["log2fc"].astype(float)
                knockout = data.get("knockout", "KO")

                # Diverging bar chart (top 20)
                df_top = df_viz.reindex(df_viz["log2fc"].abs().sort_values(ascending=False).index[:20])
                df_top = df_top.sort_values("log2fc")
                fig = go.Figure(go.Bar(
                    x=df_top["log2fc"], y=df_top["affected_gene"], orientation="h",
                    marker_color=["#3b82f6" if fc > 0 else "#f59e0b" for fc in df_top["log2fc"]],
                    text=[f"{fc:.2f}" for fc in df_top["log2fc"]], textposition="outside",
                ))
                fig.update_layout(
                    title=f"Top DE Genes \u2014 {knockout} Knockout (Mouse Brain)",
                    xaxis_title="log2 Fold Change",
                    height=max(350, len(df_top) * 25),
                    margin=dict(l=120, r=50, t=50, b=40),
                    template="plotly_dark",
                )
                st.plotly_chart(fig, use_container_width=True)

                # Volcano plot (log2FC vs -log10 pval)
                if "pvals_adj" in df_viz.columns:
                    df_viz["pvals_adj"] = df_viz["pvals_adj"].astype(float)
                    df_vol = df_viz[df_viz["pvals_adj"] > 0].copy()
                    df_vol["-log10p"] = -np.log10(df_vol["pvals_adj"].clip(lower=1e-300))
                    df_vol["significant"] = (df_vol["pvals_adj"] < 0.05) & (df_vol["log2fc"].abs() > 0.5)
                    fig_vol = go.Figure()
                    # Non-significant
                    ns = df_vol[~df_vol["significant"]]
                    fig_vol.add_trace(go.Scatter(
                        x=ns["log2fc"], y=ns["-log10p"], mode="markers",
                        marker=dict(size=5, color="gray", opacity=0.4), name="NS",
                    ))
                    # Significant
                    sig = df_vol[df_vol["significant"]]
                    fig_vol.add_trace(go.Scatter(
                        x=sig["log2fc"], y=sig["-log10p"], mode="markers",
                        marker=dict(size=7, color=["#3b82f6" if fc > 0 else "#f59e0b" for fc in sig["log2fc"]]),
                        text=sig["affected_gene"], name="Significant",
                    ))
                    fig_vol.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="red", opacity=0.5)
                    fig_vol.add_vline(x=0.5, line_dash="dash", line_color="gray", opacity=0.3)
                    fig_vol.add_vline(x=-0.5, line_dash="dash", line_color="gray", opacity=0.3)
                    fig_vol.update_layout(
                        title=f"Volcano Plot \u2014 {knockout} Knockout",
                        xaxis_title="log2 Fold Change", yaxis_title="-log10(adj p-value)",
                        height=450, template="plotly_dark",
                    )
                    st.plotly_chart(fig_vol, use_container_width=True)

                # Dot plot (cell type x gene) if multiple cell types
                if "cell_type" in df_viz.columns and df_viz["cell_type"].nunique() > 1:
                    top_genes_dot = df_viz.groupby("affected_gene")["log2fc"].apply(lambda x: x.abs().max()).nlargest(10).index.tolist()
                    df_dot = df_viz[df_viz["affected_gene"].isin(top_genes_dot)].copy()
                    df_dot["abs_fc"] = df_dot["log2fc"].abs()
                    fig_dot = go.Figure(go.Scatter(
                        x=df_dot["affected_gene"], y=df_dot["cell_type"],
                        mode="markers",
                        marker=dict(
                            size=df_dot["abs_fc"] * 8,
                            color=df_dot["log2fc"],
                            colorscale="RdBu_r", cmid=0,
                            showscale=True, colorbar=dict(title="log2FC"),
                            line=dict(width=0.5, color="white"),
                        ),
                    ))
                    fig_dot.update_layout(
                        title=f"{knockout} KO \u2014 Effects by Cell Type x Gene",
                        xaxis_title="Affected Gene", yaxis_title="Cell Type",
                        height=max(400, df_dot["cell_type"].nunique() * 30),
                        template="plotly_dark",
                    )
                    st.plotly_chart(fig_dot, use_container_width=True)

        # ══════════════════════════════════════════════════════════════════
        # 5. CROSS-SPECIES COMPARISON - Bar + Parallel Coordinates
        # ══════════════════════════════════════════════════════════════════
        elif tool == "compare_cross_species":
            mouse_data = data.get("mouse", {}).get("results", [])
            human_data = data.get("human", {}).get("results", [])

            if mouse_data and len(mouse_data) > 2:
                df_m = pd.DataFrame(mouse_data[:15])
                if "n_sig_effects" in df_m.columns:
                    df_m["n_sig_effects"] = df_m["n_sig_effects"].astype(int)
                    df_m = df_m.sort_values("n_sig_effects", ascending=True)
                    fig = go.Figure(go.Bar(
                        x=df_m["n_sig_effects"], y=df_m["cell_type"], orientation="h",
                        marker_color="#06b6d4",
                    ))
                    fig.update_layout(
                        title=f"{data.get('gene', '')} KO \u2014 Effects by Cell Type (Mouse)",
                        xaxis_title="# Significant DE Genes",
                        height=max(300, len(df_m) * 28),
                        margin=dict(l=200, r=50, t=50, b=40),
                        template="plotly_dark",
                    )
                    st.plotly_chart(fig, use_container_width=True)

            # Human CRISPRbrain results as bar
            if human_data and len(human_data) > 1:
                df_h = pd.DataFrame(human_data[:15])
                if "phenotype_score" in df_h.columns:
                    df_h["phenotype_score"] = df_h["phenotype_score"].astype(float)
                    df_h["label"] = df_h.apply(lambda r: f"{r.get('screen_name', '')[:25]} ({r.get('phenotype_name', '')})", axis=1)
                    df_h = df_h.sort_values("phenotype_score")
                    fig_h = go.Figure(go.Bar(
                        x=df_h["phenotype_score"], y=df_h["label"], orientation="h",
                        marker_color=["#ef4444" if s < 0 else "#22c55e" for s in df_h["phenotype_score"]],
                    ))
                    fig_h.update_layout(
                        title=f"{data.get('gene', '')} \u2014 Human CRISPRbrain Phenotype Scores",
                        xaxis_title="Phenotype Score",
                        height=max(300, len(df_h) * 28),
                        margin=dict(l=250, r=50, t=50, b=40),
                        template="plotly_dark",
                    )
                    st.plotly_chart(fig_h, use_container_width=True)

        # ══════════════════════════════════════════════════════════════════
        # 6. FIND KNOCKOUTS AFFECTING GENE - Network-style bar
        # ══════════════════════════════════════════════════════════════════
        elif tool == "find_knockouts_affecting_gene" and len(results) > 2:
            df_viz = pd.DataFrame(results)
            if "log2fc" in df_viz.columns and "knockout" in df_viz.columns:
                df_viz["log2fc"] = df_viz["log2fc"].astype(float)
                target_gene = data.get("target_gene", "gene")
                df_top = df_viz.reindex(df_viz["log2fc"].abs().sort_values(ascending=False).index[:20])
                df_top = df_top.sort_values("log2fc")
                fig = go.Figure(go.Bar(
                    x=df_top["log2fc"], y=df_top["knockout"], orientation="h",
                    marker_color=["#10b981" if fc > 0 else "#f43f5e" for fc in df_top["log2fc"]],
                    text=[f"{fc:.2f}" for fc in df_top["log2fc"]], textposition="outside",
                ))
                fig.update_layout(
                    title=f"Upstream Regulators of {target_gene} (KOs that affect it)",
                    xaxis_title=f"Effect on {target_gene} (log2FC)",
                    height=max(350, len(df_top) * 25),
                    margin=dict(l=120, r=50, t=50, b=40),
                    template="plotly_dark",
                )
                st.plotly_chart(fig, use_container_width=True)

        # ════════════════════════════════════════════════════════════════════
        # 7. gnomAD CONSTRAINT - pLI/LOEUF horizontal bar
        # ════════════════════════════════════════════════════════════════════
        elif tool == "query_gnomad_constraint" and len(results) >= 1:
            genes = [r.get("gene_symbol", "") for r in results]
            pli_vals = [float(r.get("pLI") or r.get("pli") or 0) for r in results]
            loeuf_vals = [float(r.get("LOEUF") or r.get("loeuf") or 0) for r in results]
            if any(v > 0 for v in pli_vals):
                fig = go.Figure()
                fig.add_trace(go.Bar(y=genes, x=pli_vals, name="pLI", orientation="h", marker_color="#8b5cf6"))
                fig.add_trace(go.Bar(y=genes, x=loeuf_vals, name="LOEUF", orientation="h", marker_color="#f59e0b"))
                fig.add_vline(x=0.9, line_dash="dash", line_color="red", opacity=0.6,
                             annotation_text="pLI=0.9 (LoF intolerant)", annotation_position="top right")
                fig.update_layout(
                    title="Gene Constraint Scores (gnomAD)",
                    xaxis_title="Score", barmode="group",
                    height=max(300, len(genes) * 50),
                    margin=dict(l=100, r=50, t=50, b=40),
                    template="plotly_dark",
                )
                st.plotly_chart(fig, use_container_width=True)

        # ════════════════════════════════════════════════════════════════════
        # 8. OpenTargets ASSOCIATIONS - Horizontal bar (score)
        # ════════════════════════════════════════════════════════════════════
        elif tool == "query_opentargets_associations" and len(results) > 2:
            import re as _re
            diseases = []
            scores = []
            for r in results[:20]:
                disease = r.get("disease", "")
                title = r.get("title", "")
                # Extract score from title like "GENE - Disease (score: 0.837)"
                m = _re.search(r"score:\s*([\d.]+)", title)
                if m and disease:
                    diseases.append(disease[:40])
                    scores.append(float(m.group(1)))
            if scores:
                # Sort by score
                paired = sorted(zip(scores, diseases), reverse=False)
                scores, diseases = zip(*paired)
                gene = data.get("gene", "Gene")
                fig = go.Figure(go.Bar(
                    x=list(scores), y=list(diseases), orientation="h",
                    marker_color=px.colors.sample_colorscale("Viridis", [s for s in scores]),
                    text=[f"{s:.3f}" for s in scores], textposition="outside",
                ))
                fig.update_layout(
                    title=f"OpenTargets Disease Associations — {gene}",
                    xaxis_title="Association Score", xaxis_range=[0, 1],
                    height=max(350, len(diseases) * 28),
                    margin=dict(l=250, r=50, t=50, b=40),
                    template="plotly_dark",
                )
                st.plotly_chart(fig, use_container_width=True)

        # ════════════════════════════════════════════════════════════════════
        # 9. Reactome/KEGG PATHWAYS - Horizontal bar (count)
        # ════════════════════════════════════════════════════════════════════
        elif tool in ("query_reactome_pathways", "query_kegg_pathways") and len(results) > 2:
            source = "Reactome" if "reactome" in tool else "KEGG"
            gene = data.get("gene", "Gene")
            pathways = [r.get("title", "").replace(f"{gene} in ", "")[:50] for r in results[:15]]
            # Simple count-based bar (each pathway = 1 membership)
            fig = go.Figure(go.Bar(
                y=pathways, x=[1] * len(pathways), orientation="h",
                marker_color="#06b6d4" if source == "Reactome" else "#10b981",
            ))
            fig.update_layout(
                title=f"{gene} — {source} Pathway Membership",
                xaxis_title="", xaxis_showticklabels=False,
                height=max(300, len(pathways) * 28),
                margin=dict(l=300, r=50, t=50, b=40),
                template="plotly_dark",
            )
            st.plotly_chart(fig, use_container_width=True)

        # ════════════════════════════════════════════════════════════════════
        # 10. CIViC EVIDENCE - Grouped by evidence level
        # ════════════════════════════════════════════════════════════════════
        elif tool == "query_civic_evidence" and len(results) > 1:
            import re as _re
            gene = data.get("gene", "Gene")
            level_counts = {}
            for r in results:
                title = r.get("title", "")
                m = _re.search(r"Level\s+(\w)", title)
                lvl = m.group(1) if m else "?"
                level_counts[lvl] = level_counts.get(lvl, 0) + 1
            if level_counts:
                levels = sorted(level_counts.keys())
                counts = [level_counts[l] for l in levels]
                colors = {"A": "#22c55e", "B": "#84cc16", "C": "#eab308", "D": "#f97316", "E": "#ef4444"}
                fig = go.Figure(go.Bar(
                    x=levels, y=counts,
                    marker_color=[colors.get(l, "#6b7280") for l in levels],
                    text=counts, textposition="outside",
                ))
                fig.update_layout(
                    title=f"{gene} — CIViC Evidence by Level (A=validated, E=preclinical)",
                    xaxis_title="Evidence Level", yaxis_title="# Evidence Items",
                    height=350, template="plotly_dark",
                )
                st.plotly_chart(fig, use_container_width=True)

        # ════════════════════════════════════════════════════════════════════
        # 11. UniProt PROTEIN FUNCTION - Domain summary
        # ════════════════════════════════════════════════════════════════════
        elif tool == "query_protein_function" and len(results) >= 1:
            gene = data.get("gene", "Gene")
            # Show a simple info card with key fields
            for r in results[:3]:
                title = r.get("title", "")
                summary = r.get("summary", "")
                if "Function:" in summary and len(summary) > 50:
                    st.info(f"**{title}**\n\n{summary[:500]}")

        # ════════════════════════════════════════════════════════════════════
        # 12. Human Protein Atlas - Tissue expression bar
        # ════════════════════════════════════════════════════════════════════
        elif tool == "query_protein_expression" and len(results) > 2:
            gene = data.get("gene", "Gene")
            tissues = []
            levels = []
            level_map = {"Not detected": 0, "Low": 1, "Medium": 2, "High": 3}
            for r in results:
                title = r.get("title", "")
                if " in " in title and ":" in title:
                    tissue = title.split(" in ")[1].split(":")[0].strip()
                    level_str = title.split(":")[-1].strip()
                    tissues.append(tissue)
                    levels.append(level_map.get(level_str, 0))
            if tissues:
                paired = sorted(zip(levels, tissues), reverse=True)
                levels_s, tissues_s = zip(*paired)
                color_map = {0: "#374151", 1: "#fbbf24", 2: "#f97316", 3: "#ef4444"}
                fig = go.Figure(go.Bar(
                    y=list(tissues_s)[:20], x=list(levels_s)[:20], orientation="h",
                    marker_color=[color_map.get(l, "#6b7280") for l in levels_s[:20]],
                ))
                fig.update_layout(
                    title=f"{gene} Protein Expression (Human Protein Atlas)",
                    xaxis_title="Expression Level",
                    xaxis=dict(tickvals=[0, 1, 2, 3], ticktext=["Not detected", "Low", "Medium", "High"]),
                    height=max(350, min(len(tissues), 20) * 28),
                    margin=dict(l=200, r=50, t=50, b=40),
                    template="plotly_dark",
                )
                st.plotly_chart(fig, use_container_width=True)

        # ════════════════════════════════════════════════════════════════════
        # 13. Monarch ASSOCIATIONS - Disease dot plot
        # ════════════════════════════════════════════════════════════════════
        elif tool == "query_monarch_associations" and len(results) > 2:
            gene = data.get("gene", "Gene")
            diseases = [r.get("disease", r.get("title", ""))[:45] for r in results if r.get("disease") or r.get("title")]
            if diseases:
                fig = go.Figure(go.Scatter(
                    x=[1] * len(diseases[:15]), y=diseases[:15],
                    mode="markers",
                    marker=dict(size=14, color="#a855f7", symbol="diamond",
                               line=dict(width=1, color="white")),
                ))
                fig.update_layout(
                    title=f"{gene} — Monarch Disease/Phenotype Associations",
                    xaxis_showticklabels=False, xaxis_title="",
                    height=max(300, len(diseases[:15]) * 30),
                    margin=dict(l=300, r=50, t=50, b=40),
                    template="plotly_dark",
                )
                st.plotly_chart(fig, use_container_width=True)

        # ════════════════════════════════════════════════════════════════════
        # 14. gnomAD VARIANTS - Bar by clinical significance
        # ════════════════════════════════════════════════════════════════════
        elif tool == "query_gnomad_variants" and len(results) > 1:
            import re as _re
            gene = data.get("gene", "Gene")
            sig_counts = {}
            for r in results:
                title = r.get("title", "")
                for sig in ["Pathogenic", "Likely pathogenic", "Uncertain significance", "Likely benign", "Benign"]:
                    if sig.lower() in title.lower():
                        sig_counts[sig] = sig_counts.get(sig, 0) + 1
                        break
                else:
                    sig_counts["Other"] = sig_counts.get("Other", 0) + 1
            if sig_counts:
                sigs = list(sig_counts.keys())
                counts = list(sig_counts.values())
                sig_colors = {"Pathogenic": "#ef4444", "Likely pathogenic": "#f97316",
                             "Uncertain significance": "#eab308", "Likely benign": "#84cc16",
                             "Benign": "#22c55e", "Other": "#6b7280"}
                fig = go.Figure(go.Bar(
                    x=sigs, y=counts,
                    marker_color=[sig_colors.get(s, "#6b7280") for s in sigs],
                    text=counts, textposition="outside",
                ))
                fig.update_layout(
                    title=f"{gene} — ClinVar Variant Classifications (gnomAD)",
                    xaxis_title="Clinical Significance", yaxis_title="# Variants",
                    height=350, template="plotly_dark",
                )
                st.plotly_chart(fig, use_container_width=True)

        # ════════════════════════════════════════════════════════════════════
        # 15. OpenTargets TRACTABILITY - Modality pie chart
        # ════════════════════════════════════════════════════════════════════
        elif tool == "query_opentargets_tractability" and len(results) >= 1:
            gene = data.get("gene", "Gene")
            modalities = []
            for r in results:
                title = r.get("title", "")
                # Extract modality from title like "GENE Tractability: Small molecule"
                if ":" in title:
                    mod = title.split(":")[-1].strip()
                    modalities.append(mod)
            if modalities:
                mod_counts = {}
                for m in modalities:
                    mod_counts[m] = mod_counts.get(m, 0) + 1
                fig = go.Figure(go.Pie(
                    labels=list(mod_counts.keys()), values=list(mod_counts.values()),
                    hole=0.4,
                    marker_colors=["#8b5cf6", "#06b6d4", "#f59e0b", "#ef4444", "#10b981", "#6b7280"],
                ))
                fig.update_layout(
                    title=f"{gene} — Druggability Modalities (OpenTargets)",
                    height=400, template="plotly_dark",
                )
                st.plotly_chart(fig, use_container_width=True)

        # ════════════════════════════════════════════════════════════════════
        # 16. cBioPortal MUTATIONS - Bar by disease/cancer type
        # ════════════════════════════════════════════════════════════════════
        elif tool == "query_cbioportal_mutations" and len(results) > 1:
            gene = data.get("gene", "Gene")
            disease_counts = {}
            for r in results:
                disease = r.get("disease", "Unknown")[:35]
                if disease:
                    disease_counts[disease] = disease_counts.get(disease, 0) + 1
            if disease_counts:
                sorted_d = sorted(disease_counts.items(), key=lambda x: x[1], reverse=True)[:15]
                diseases, counts = zip(*sorted_d)
                fig = go.Figure(go.Bar(
                    y=list(diseases), x=list(counts), orientation="h",
                    marker_color="#f43f5e",
                    text=list(counts), textposition="outside",
                ))
                fig.update_layout(
                    title=f"{gene} — Cancer Mutation Frequency (cBioPortal)",
                    xaxis_title="# Mutations", yaxis=dict(autorange="reversed"),
                    height=max(300, len(diseases) * 30),
                    margin=dict(l=250, r=50, t=50, b=40),
                    template="plotly_dark",
                )
                st.plotly_chart(fig, use_container_width=True)

        # ════════════════════════════════════════════════════════════════════
        # 17. NCBI GTR - Genetic test count bar
        # ════════════════════════════════════════════════════════════════════
        elif tool == "query_genetic_tests" and len(results) > 1:
            gene = data.get("gene", "Gene")
            test_types = {}
            for r in results:
                title = r.get("title", "")
                # Categorize by test type keywords
                if "sequencing" in title.lower():
                    test_types["Sequencing"] = test_types.get("Sequencing", 0) + 1
                elif "panel" in title.lower():
                    test_types["Gene Panel"] = test_types.get("Gene Panel", 0) + 1
                elif "deletion" in title.lower() or "dup" in title.lower():
                    test_types["Del/Dup"] = test_types.get("Del/Dup", 0) + 1
                else:
                    test_types["Other"] = test_types.get("Other", 0) + 1
            if test_types:
                fig = go.Figure(go.Bar(
                    x=list(test_types.keys()), y=list(test_types.values()),
                    marker_color="#06b6d4",
                    text=list(test_types.values()), textposition="outside",
                ))
                fig.update_layout(
                    title=f"{gene} — Available Genetic Tests (NCBI GTR)",
                    xaxis_title="Test Type", yaxis_title="# Tests",
                    height=350, template="plotly_dark",
                )
                st.plotly_chart(fig, use_container_width=True)

        # ════════════════════════════════════════════════════════════════════
        # 18. Human CRISPRbrain SCREEN HITS - Phenotype score bar
        # ════════════════════════════════════════════════════════════════════
        elif tool == "query_human_screen_hits" and len(results) > 1:
            df_viz = pd.DataFrame(results)
            if "phenotype_score" in df_viz.columns:
                df_viz["phenotype_score"] = df_viz["phenotype_score"].astype(float)
                gene = data.get("gene", "Gene")
                df_viz["label"] = df_viz.apply(
                    lambda r: f"{str(r.get('screen_name', ''))[:25]} ({r.get('phenotype_name', '')})", axis=1
                )
                df_viz = df_viz.sort_values("phenotype_score")
                fig = go.Figure(go.Bar(
                    x=df_viz["phenotype_score"], y=df_viz["label"], orientation="h",
                    marker_color=["#ef4444" if s < 0 else "#22c55e" for s in df_viz["phenotype_score"]],
                    text=[f"{s:.2f}" for s in df_viz["phenotype_score"]], textposition="outside",
                ))
                fig.update_layout(
                    title=f"{gene} — CRISPRbrain Phenotype Scores",
                    xaxis_title="Phenotype Score (negative = depleted, positive = enriched)",
                    height=max(300, len(df_viz) * 28),
                    margin=dict(l=250, r=50, t=50, b=40),
                    template="plotly_dark",
                )
                st.plotly_chart(fig, use_container_width=True)

        # ════════════════════════════════════════════════════════════════════
        # 19. Human Expression (Cohorts) - Tissue bar
        # ════════════════════════════════════════════════════════════════════
        elif tool == "query_human_expression" and len(results) > 2:
            df_viz = pd.DataFrame(results)
            gene = data.get("gene", "Gene")
            if "mean_expr" in df_viz.columns and "tissue" in df_viz.columns:
                df_viz["mean_expr"] = df_viz["mean_expr"].astype(float)
                df_viz = df_viz.sort_values("mean_expr", ascending=True)
                fig = go.Figure(go.Bar(
                    x=df_viz["mean_expr"], y=df_viz["tissue"], orientation="h",
                    marker_color="#a855f7",
                    text=[f"{v:.2f}" for v in df_viz["mean_expr"]], textposition="outside",
                ))
                fig.update_layout(
                    title=f"{gene} Expression (BioFINDER/ROSMAP Cohorts)",
                    xaxis_title="Mean Expression",
                    height=max(300, len(df_viz) * 28),
                    margin=dict(l=200, r=50, t=50, b=40),
                    template="plotly_dark",
                )
                st.plotly_chart(fig, use_container_width=True)

        # ════════════════════════════════════════════════════════════════════
        # 20. Tahoe-100M DRUG CLUSTERS - Cluster distribution bar
        # ════════════════════════════════════════════════════════════════════
        elif tool == "query_tahoe_drug_clusters" and len(results) > 2:
            df_viz = pd.DataFrame(results)
            drug = data.get("drug", "Drug")
            if "cluster" in df_viz.columns and "cell_count" in df_viz.columns:
                df_viz["cell_count"] = df_viz["cell_count"].astype(int)
                df_viz = df_viz.sort_values("cell_count", ascending=False)[:20]
                # Color by cell line if present
                if "cell_line" in df_viz.columns and df_viz["cell_line"].nunique() > 1:
                    fig = px.bar(
                        df_viz, x="cluster", y="cell_count", color="cell_line",
                        barmode="group",
                    )
                else:
                    fig = go.Figure(go.Bar(
                        x=df_viz["cluster"], y=df_viz["cell_count"],
                        marker_color="#06b6d4",
                    ))
                fig.update_layout(
                    title=f"{drug} — Transcriptomic Cluster Distribution (Tahoe-100M)",
                    xaxis_title="Cluster", yaxis_title="Cell Count",
                    height=400, template="plotly_dark",
                )
                st.plotly_chart(fig, use_container_width=True)

        # ════════════════════════════════════════════════════════════════════
        # 21. Tahoe-100M MOA DISTRIBUTION - Horizontal bar
        # ════════════════════════════════════════════════════════════════════
        elif tool == "query_tahoe_moa_distribution" and len(results) > 2:
            df_viz = pd.DataFrame(results)
            if "moa" in df_viz.columns and "drug_count" in df_viz.columns:
                df_viz["drug_count"] = df_viz["drug_count"].astype(int)
                df_viz = df_viz.sort_values("drug_count", ascending=True).tail(20)
                fig = go.Figure(go.Bar(
                    x=df_viz["drug_count"], y=df_viz["moa"], orientation="h",
                    marker_color="#10b981",
                    text=df_viz["drug_count"], textposition="outside",
                ))
                fig.update_layout(
                    title="Tahoe-100M — Compounds by Mechanism of Action",
                    xaxis_title="# Compounds",
                    height=max(350, len(df_viz) * 28),
                    margin=dict(l=250, r=50, t=50, b=40),
                    template="plotly_dark",
                )
                st.plotly_chart(fig, use_container_width=True)

        # ════════════════════════════════════════════════════════════════════
        # 22. NDD RISK GENES - FDR score bar
        # ════════════════════════════════════════════════════════════════════
        elif tool == "list_ndd_risk_genes" and len(results) > 3:
            df_viz = pd.DataFrame(results)
            if "fdr" in df_viz.columns and "gene" in df_viz.columns:
                df_viz["fdr"] = df_viz["fdr"].astype(float)
                df_viz["-log10fdr"] = -np.log10(df_viz["fdr"].clip(lower=1e-300))
                df_viz = df_viz.sort_values("-log10fdr", ascending=True).tail(25)
                fig = go.Figure(go.Bar(
                    x=df_viz["-log10fdr"], y=df_viz["gene"], orientation="h",
                    marker_color="#f97316",
                    text=[f"{f:.1e}" for f in df_viz["fdr"]], textposition="outside",
                ))
                fig.add_vline(x=-np.log10(0.05), line_dash="dash", line_color="red",
                             annotation_text="FDR=0.05", annotation_position="top right")
                fig.update_layout(
                    title="NDD/ASD Risk Genes (TADA analysis)",
                    xaxis_title="-log10(FDR)",
                    height=max(350, len(df_viz) * 25),
                    margin=dict(l=100, r=80, t=50, b=40),
                    template="plotly_dark",
                )
                st.plotly_chart(fig, use_container_width=True)

        # ════════════════════════════════════════════════════════════════════
        # 23. BRAIN CELL TYPES - Cell count horizontal bar
        # ════════════════════════════════════════════════════════════════════
        elif tool == "list_brain_cell_types" and len(results) > 3:
            df_viz = pd.DataFrame(results)
            if "cell_count" in df_viz.columns and "cell_type" in df_viz.columns:
                df_viz["cell_count"] = df_viz["cell_count"].astype(int)
                df_viz = df_viz.sort_values("cell_count", ascending=True).tail(25)
                fig = go.Figure(go.Bar(
                    x=df_viz["cell_count"], y=df_viz["cell_type"], orientation="h",
                    marker_color="#06b6d4",
                ))
                fig.update_layout(
                    title="Mouse Brain Cell Types (CRISPR Atlas)",
                    xaxis_title="# Cells",
                    height=max(350, len(df_viz) * 25),
                    margin=dict(l=200, r=50, t=50, b=40),
                    template="plotly_dark",
                )
                st.plotly_chart(fig, use_container_width=True)

    # ══════════════════════════════════════════════════════════════════════
    # MULTI-GENE HEATMAP (if multiple expression queries were made)
    # ══════════════════════════════════════════════════════════════════════
    if len(expression_results) > 1:
        # Build matrix: genes x brain regions
        all_tissues = expression_results[0]["tissues"]
        genes = [er["gene"] for er in expression_results]
        matrix = []
        for er in expression_results:
            tissue_map = dict(zip(er["tissues"], er["tpms"]))
            matrix.append([tissue_map.get(t, 0) for t in all_tissues])

        fig_heat = go.Figure(go.Heatmap(
            z=matrix, x=all_tissues, y=genes,
            colorscale="Viridis",
            text=np.round(matrix, 2), texttemplate="%{text}",
        ))
        fig_heat.update_layout(
            title="Multi-Gene Expression Heatmap (Brain Regions)",
            height=max(300, len(genes) * 40 + 100),
            margin=dict(l=100, r=50, t=50, b=120),
            xaxis=dict(tickangle=45),
            template="plotly_dark",
        )
        st.plotly_chart(fig_heat, use_container_width=True)


# ── Load icon once for reuse ──
import base64 as _b64
_icon_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "neuroplex_icon.svg")
with open(_icon_path, "r") as _f:
    _icon_b64 = _b64.b64encode(_f.read().encode()).decode()

# ── Streamlit UI ──
with st.sidebar:
    st.markdown(
        f'<div style="display:flex;align-items:center;gap:10px;margin-bottom:0.25rem;">'
        f'<img src="data:image/svg+xml;base64,{_icon_b64}" width="32" height="32" style="border-radius:6px;">'
        f'<span style="font-size:1.5rem;font-weight:700;">NeuroPlex</span>'
        f'</div>',
        unsafe_allow_html=True,
    )
    st.caption(_t("subtitle"))
    
    # ── Language Toggle ──
    if "lang" not in st.session_state:
        st.session_state["lang"] = "en"
    _lang_options = {"en": "\U0001F1FA\U0001F1F8 English", "ja": "\U0001F1EF\U0001F1F5 \u65e5\u672c\u8a9e"}
    _selected_lang = st.radio(
        "Language / \u8a00\u8a9e",
        options=list(_lang_options.keys()),
        format_func=lambda k: _lang_options[k],
        index=list(_lang_options.keys()).index(st.session_state["lang"]),
        horizontal=True,
        key="_lang_radio",
    )
    if _selected_lang != st.session_state["lang"]:
        st.session_state["lang"] = _selected_lang
        st.rerun()
    
    if _sql_ok:
        st.success(_t("sql_ready", msg=_sql_msg), icon="\u2705")
    else:
        st.error(_t("sql_error", msg=_sql_msg), icon="\u26A0\uFE0F")
    
    st.divider()
    
    # ── Dataset Focus ──
    st.subheader(_t("dataset_focus"))

    _DATASET_OPTS = [
        ("mouse",      _t("ds_mouse")),
        ("human",      _t("ds_human")),
        ("expression", _t("ds_expression")),
        ("pharma",     _t("ds_pharma")),
        ("tahoe",      _t("ds_tahoe")),
        ("genetics",   _t("ds_genetics")),
        ("druggability", _t("ds_druggability")),
        ("protein",    _t("ds_protein")),
        ("pathways",   _t("ds_pathways")),
    ]
    _ALL_DS_KEYS = [k for k, _ in _DATASET_OPTS]

    # Initialize all datasets as selected by default
    for _k, _ in _DATASET_OPTS:
        if f"ds_{_k}" not in st.session_state:
            st.session_state[f"ds_{_k}"] = True

    _dc1, _dc2 = st.columns(2)
    with _dc1:
        if st.button(_t("btn_all"), key="ds_btn_all", use_container_width=True):
            for _k, _ in _DATASET_OPTS:
                st.session_state[f"ds_{_k}"] = True
            st.rerun()
    with _dc2:
        if st.button(_t("btn_none"), key="ds_btn_none", use_container_width=True):
            for _k, _ in _DATASET_OPTS:
                st.session_state[f"ds_{_k}"] = False
            st.rerun()

    active_datasets = []
    for _key, _label in _DATASET_OPTS:
        if st.checkbox(_label, key=f"ds_{_key}"):
            active_datasets.append(_key)
    if not active_datasets:
        active_datasets = _ALL_DS_KEYS

    # Dataset info as collapsible expanders
    _DATASET_INFO = {
        "mouse":      ("\U0001F42D Mouse (CRISPR Atlas + IMPC)",
                       "**PerturbAI CRISPR Atlas**\n- 7.7M cells, in vivo mouse brain\n- 2,046 gene knockouts\n- 745M DE results\n- Hundreds of neuron types\n\n**IMPC**\n- Whole-organism knockout phenotyping\n- MP ontology (locomotor, cardiac, metabolic, neurological)\n- Orthology: mouse\u2194human gene mappings\n- Allele registry: live lines + ordering links\n- Publications: semantic search"),
        "human":      ("\U0001F9EC Human CRISPRbrain",
                       "- iPSC neurons, microglia, astrocytes\n- Genome-wide CRISPRi/a\n- Survival, tau, TDP-43, phagocytosis\n- MAPT-V337M, PSAP KO models"),
        "expression": ("\U0001F4CA Expression (GTEx + Cohorts)",
                       "- GTEx v8: 13 brain subregions (baseline TPM)\n- BioFINDER + ROSMAP cohorts\n- HCRT, HCRTR1, HCRTR2 localization"),
        "pharma":     ("\U0001F48A Pharmacology (LINCS + ChEMBL)",
                       "**LINCS L1000**\n- Orexin antagonist signatures\n- SH-SY5Y, NPC neuronal lines\n- Suvorexant, lemborexant, daridorexant\n\n"
                       "**ChEMBL**\n- OX1R / OX2R binding (Ki, IC50)\n- Clinical DORAs + tool compounds\n- Selectivity profiles"),
        "tahoe":      ("\U0001F52C Tahoe-100M Drug Atlas",
                       "- ~100M cells, drug perturbation scRNA-seq\n- Diverse compounds across cell lines\n- MOA annotations + cluster assignments\n- SMILES + PubChem CIDs"),
        "genetics":   ("\U0001F9EC Genetics (gnomAD + CIViC + GTR)",
                       "**gnomAD v4**\n- Gene constraint: pLI, LOEUF, missense Z\n- Loss-of-function intolerance scores\n\n"
                       "**CIViC**\n- Clinical variant evidence (Level A\u2013E)\n- Therapy associations & significance\n\n"
                       "**cBioPortal**\n- Cancer somatic mutations\n- Protein changes & study context\n\n"
                       "**NCBI GTR**\n- Available clinical genetic tests"),
        "druggability": ("\U0001F3AF Druggability (OpenTargets)",
                       "- Target-disease association scores\n- Evidence from genetics, literature, pathways\n- Tractability: small molecule, antibody, PROTAC\n- Therapeutic area classification"),
        "protein":    ("\U0001F52C Protein (HPA + Monarch)",
                       "**Human Protein Atlas**\n- Tissue-level protein expression\n- Subcellular localization\n- Brain, liver, kidney, etc.\n\n"
                       "**Monarch Initiative**\n- Gene-disease associations\n- Phenotype matching from model organisms\n- Cross-species phenotypes"),
        "pathways":   ("\U0001F5FA\uFE0F Pathways (UniProt + Reactome + KEGG)",
                       "**UniProt**\n- Protein function & domain architecture\n- Subcellular location & disease links\n\n"
                       "**Reactome**\n- Biological pathway hierarchy\n- Reactions & compartments\n\n"
                       "**KEGG**\n- Metabolic & signaling pathway maps\n- Pathway classes & interactions"),
    }
    with st.expander("📋 Dataset Summaries", expanded=False):
        for _ds_key in active_datasets:
            if _ds_key in _DATASET_INFO:
                _title, _body = _DATASET_INFO[_ds_key]
                st.markdown(f"**{_title}**")
                st.markdown(_body)
                st.markdown("---")

    st.divider()

    # ── Documentation ──
    st.subheader("📚 Documentation")
    _LITERATURE_DOCS = [
        ("PubMed", "Article search, PubTator annotations, and PMC full-text handoff", "https://pubmed.ncbi.nlm.nih.gov/"),
        ("ClinicalTrials.gov", "Recruiting-study search, eligibility text, and site details", "https://clinicaltrials.gov/"),
        ("ClinVar", "Clinical significance and review-status context for variants", "https://www.ncbi.nlm.nih.gov/clinvar/"),
        ("OpenFDA", "FAERS, recalls, device events, labels, and U.S. approval context", "https://open.fda.gov/"),
        ("Semantic Scholar", "TLDRs, citation graphs, references, and recommendations", "https://www.semanticscholar.org/"),
        ("EMA", "EU regulatory, safety, and shortage context for medicines", "https://www.ema.europa.eu/"),
        ("MedlinePlus", "Plain-language disease/symptom context for discover and disease clinical features", "https://medlineplus.gov/"),
        ("PharmGKB / CPIC", "Pharmacogenomic recommendations, frequencies, and clinical annotations", "https://www.pharmgkb.org/"),
    ]
    for _name, _desc, _url in _LITERATURE_DOCS:
        with st.expander(_name, expanded=False):
            st.markdown(f"{_desc}\n\n[{_name}]({_url})")

    st.divider()
    if st.button("\U0001F5D1\uFE0F  Clear Chat", use_container_width=True):
        st.session_state.messages = []
        st.rerun()
    
    st.markdown("---")
    st.markdown("**Example queries:**")
    
    _focus_set = set(active_datasets)
    if _focus_set == {"mouse"}:
        examples = [
            "What happens when Psen1 is knocked out?",
            "Which KOs affect Hcrt expression in neurons?",
            "What IMPC phenotypes does Hcrt knockout show?",
            "Find the human ortholog of Psen1 via IMPC",
            "Are there live Hcrt mouse lines available to order?",
        ]
    elif _focus_set == {"human"}:
        examples = [
            "What screens show PSEN1 as a hit?",
            "List human neuron screens",
            "Is SOD2 protective for neuron survival?",
            "What is MAPT expression in brain?",
        ]
    elif _focus_set == {"tahoe"}:
        examples = [
            "What drugs are in Tahoe-100M?",
            "Show suvorexant cluster distribution",
            "List orexin antagonists in Tahoe-100M",
            "What MOA classes are represented?",
        ]
    elif _focus_set == {"pharma"}:
        examples = [
            "Compare suvorexant vs lemborexant binding",
            "Show OX2R selectivity for DORAs",
            "What genes does lemborexant perturb?",
            "Compare DORA potencies across OX1R and OX2R",
        ]
    elif _focus_set == {"genetics"}:
        examples = _t("examples_default") if st.session_state.get("lang") == "ja" else [
            "What is the constraint score for PSEN1?",
            "Show pathogenic variants in HCRTR1",
            "Is SOD1 loss-of-function intolerant?",
            "Compare pLI scores for orexin genes",
        ]
    elif _focus_set == {"druggability"}:
        examples = _t("examples_default") if st.session_state.get("lang") == "ja" else [
            "What diseases are associated with HCRT?",
            "Is HCRTR2 druggable?",
            "Show OpenTargets associations for PSEN1",
            "What is the tractability of APP?",
        ]
    else:
        examples = [
            "Compare PSEN1 across mouse and human",
            "Summarize all datasets",
            "List ASD risk genes",
            "Is APP important for neuron survival?",
            "What drugs target orexin receptors?",
        ]
    
    for ex in examples:
        if st.button(ex, key=f"ex_{ex}", use_container_width=True):
            st.session_state.messages = [{"role": "user", "content": ex}]
            st.rerun()

# ── Get user identity ──
if "user_email" not in st.session_state:
    st.session_state.user_email = _get_user_email()

# ── Main Content: Two-column layout ──
col_chat, col_history = st.columns([3, 1], gap="large")

with col_chat:
    import base64 as _b64
    _icon_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "neuroplex_icon.svg")
    with open(_icon_path, "r") as _f:
        _icon_b64 = _b64.b64encode(_f.read().encode()).decode()
    st.markdown(
        f'<div style="display:flex;align-items:center;gap:14px;margin-bottom:0.25rem;">'
        f'<img src="data:image/svg+xml;base64,{_icon_b64}" width="52" height="52" style="border-radius:8px;">'
        f'<h1 style="margin:0;padding:0;line-height:1.2;">NeuroPlex</h1>'
        f'</div>',
        unsafe_allow_html=True,
    )
    _DS_SHORT = {"mouse": _t("ds_short_mouse"), "human": _t("ds_short_human"), "expression": _t("ds_short_expression"), "pharma": _t("ds_short_pharma"), "tahoe": _t("ds_short_tahoe"), "genetics": _t("ds_short_genetics"), "druggability": _t("ds_short_druggability"), "protein": _t("ds_short_protein"), "pathways": _t("ds_short_pathways")}
    _active_label = ", ".join(_DS_SHORT.get(k, k) for k in active_datasets) if len(active_datasets) < 9 else _t("all_datasets")
    st.caption(_t("focus_caption", label=_active_label))

    if "messages" not in st.session_state:
        st.session_state.messages = []

    def _render_md(text: str):
        """Render markdown, stripping code fences the LLM may wrap around the response."""
        _stripped = text.strip()
        if _stripped.startswith("```") and _stripped.endswith("```"):
            # Remove opening fence (```markdown or ```)
            first_newline = _stripped.find("\n")
            if first_newline != -1:
                _stripped = _stripped[first_newline + 1:]
            # Remove closing fence
            _stripped = _stripped.rsplit("```", 1)[0].rstrip()
        st.markdown(_stripped, unsafe_allow_html=True)

    for msg in st.session_state.messages:
        with st.chat_message(msg["role"]):
            _render_md(msg["content"])
            # Re-render charts for assistant messages that have stored tool results
            if msg["role"] == "assistant" and msg.get("tool_results"):
                render_charts(msg["tool_results"])

    if prompt := st.chat_input(_t("chat_placeholder")):
        st.session_state.messages.append({"role": "user", "content": prompt})
        with st.chat_message("user"):
            st.markdown(prompt)

    if st.session_state.messages and st.session_state.messages[-1]["role"] == "user":
        with st.chat_message("assistant"):
            _ds_spinner = ", ".join(_DS_SHORT.get(k, k) for k in active_datasets[:3])
            with st.spinner(f"Querying {_ds_spinner}..."):
                try:
                    response = run_agent(st.session_state.messages, active_datasets=active_datasets)
                except Exception as e:
                    logger.error(f"Agent error: {traceback.format_exc()}")
                    response = f"Error: {type(e).__name__}: {e}"
            _render_md(response)
            # Render charts from tool results
            if st.session_state.get("tool_results"):
                render_charts(st.session_state.tool_results)
        st.session_state.messages.append({"role": "assistant", "content": response})
        # Log query to persistent history
        log_query(
            user_email=st.session_state.user_email,
            species_focus=",".join(active_datasets),
            query_text=st.session_state.messages[-2]["content"],
            response_text=response,
        )

# ── Right Panel: Query History ──
with col_history:
    st.markdown("### \U0001F4DC Query Log")
    st.caption("Recent queries from all users")

    _ref_col, _clr_col = st.columns(2)
    with _ref_col:
        if st.button("\U0001F504 Refresh", key="refresh_history", use_container_width=True):
            st.cache_data.clear()
            st.rerun()
    with _clr_col:
        if st.button("\U0001F5D1\uFE0F Clear Log", key="clear_log_btn", use_container_width=True, type="secondary"):
            st.session_state["confirm_clear_log"] = True

    # Confirmation dialog for clearing query log
    if st.session_state.get("confirm_clear_log"):
        st.warning("\u26A0\uFE0F **This action is irreversible.** All query history will be permanently deleted and cannot be recovered.")
        _y_col, _n_col = st.columns(2)
        with _y_col:
            if st.button("Yes, clear log", key="confirm_clear_yes", use_container_width=True, type="primary"):
                try:
                    _execute_query(f"DELETE FROM {TABLE_QUERY_LOG}", timeout=30)
                    st.cache_data.clear()
                    st.session_state["confirm_clear_log"] = False
                    st.success("Query log cleared.")
                    st.rerun()
                except Exception as _e:
                    st.error(f"Failed to clear log: {_e}")
        with _n_col:
            if st.button("Cancel", key="confirm_clear_no", use_container_width=True):
                st.session_state["confirm_clear_log"] = False
                st.rerun()

    history = fetch_query_history(limit=50)

    if not history:
        st.info("No queries logged yet. Ask a question to get started!")
    else:
        for i, entry in enumerate(history):
            ts = entry.get("timestamp", "")[:16]  # trim seconds
            user = entry.get("user_email", "unknown").split("@")[0]  # show first part
            query = entry.get("query_text", "")
            focus = entry.get("species_focus", "")
            response = entry.get("response_text", "")

            # Focus badge
            _first_focus = focus.split(",")[0] if focus else "both"
            focus_icon = {"both": "\U0001F310", "mouse": "\U0001F42D", "human": "\U0001F9EC", "pharma": "\U0001F48A", "expression": "\U0001F4CA", "tahoe": "\U0001F52C"}.get(_first_focus, "\U0001F50D")

            with st.expander(f"{focus_icon} {query[:60]}{'...' if len(query) > 60 else ''}", expanded=False):
                st.markdown(f"**{user}** \u00b7 {ts} \u00b7 {focus_icon} {focus}")
                if st.button("\u21a9 Re-run", key=f"rerun_{i}", use_container_width=True):
                    st.session_state.messages = [{"role": "user", "content": query}]
                    st.rerun()
                st.markdown("---")
                if len(query) > 60:
                    st.caption(f"\U0001F4CB {query}")
                    st.markdown("---")
                _hist_text = response[:2000] + ("..." if len(response) > 2000 else "")
                st.markdown(_hist_text, unsafe_allow_html=True)

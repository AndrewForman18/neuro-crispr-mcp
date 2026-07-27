# Databricks notebook source
# /// script
# [tool.databricks.environment]
# environment_version = "5"
# ///
# DBTITLE 1,Set VAL environment (must run first)
import os, sys, importlib
os.environ["NEUROPLEX_ENV"] = "val"

# Force reload so edits to neuroplex_env.yml are picked up in the same session
sys.path.insert(0, "/Workspace/Users/andrew_forman@eisai.com/neuro-crispr-mcp")
for mod in [k for k in sys.modules if k.startswith("config")]:
    del sys.modules[mod]
from config.neuroplex_config import load_config
CFG = load_config()
assert CFG.environment == "val", f"Expected val, got {CFG.environment}"
assert CFG.catalog == "dhbl_discovery_us_val", f"Wrong catalog: {CFG.catalog}"
print(f"✅ Environment: {CFG.environment}")
print(f"   Catalog: {CFG.catalog}.{CFG.schema}")
print(f"   Warehouse: {CFG.sql_warehouse_id}")

# COMMAND ----------

# DBTITLE 1,Step 1: CRISPRbrain screens — status
print("Running 01_ingest_crisprbrain...")

# COMMAND ----------

# DBTITLE 1,Install openpyxl (required by 01_ingest_crisprbrain)
# MAGIC %pip install openpyxl --quiet

# COMMAND ----------

# DBTITLE 1,Step 1: %run 01_ingest_crisprbrain
# MAGIC %run ./01_ingest_crisprbrain

# COMMAND ----------

# DBTITLE 1,Step 2: GTEx + LINCS + ChEMBL — status
print("Running 02_ingest_pharmacology...")

# COMMAND ----------

# DBTITLE 1,Step 2b: gnomAD ETL — status
print("Running 03_ingest_gnomad (gnomAD v4 constraint + ClinVar variants)...")
print("  Gene panel: 31 neurodegeneration / orexin genes")
print("  API: https://gnomad.broadinstitute.org/api (GraphQL, public, no auth)")
print("  Target: neuroplex_gnomad (VARIANT payload)")

# COMMAND ----------

# DBTITLE 1,Step 2: %run 02_ingest_pharmacology
# MAGIC %run ./02_ingest_pharmacology

# COMMAND ----------

# DBTITLE 1,Step 3: Create neuroplex_query_log
spark.sql(f"""
CREATE TABLE IF NOT EXISTS {CFG.catalog}.{CFG.schema}.neuroplex_query_log (
    log_id        STRING NOT NULL COMMENT 'UUID',
    ts            TIMESTAMP NOT NULL COMMENT 'Query timestamp',
    user_session  STRING COMMENT 'Session identifier',
    user_query    STRING COMMENT 'Raw user query text',
    datasets_used ARRAY<STRING> COMMENT 'Dataset keys queried',
    tool_calls    ARRAY<STRING> COMMENT 'MCP tools invoked',
    response_ms   INT COMMENT 'Response time in ms',
    environment   STRING COMMENT 'neuroplex environment'
)
USING DELTA
COMMENT 'NeuroPlex query audit log'
""")
print(f"✅ {CFG.catalog}.{CFG.schema}.neuroplex_query_log created")

# COMMAND ----------

# DBTITLE 1,Step 4: Placeholder tables for large datasets
# Schemas derived from SELECT statements in neuro_mcp_server/server.py
_PLACEHOLDER_COMMENT = "PLACEHOLDER — correct schema, no data loaded. Full ingestion pending."

placeholder_ddls = {
    "wholebrain_crispr_atlas": """
        cell_barcode        STRING  COMMENT 'Unique cell identifier',
        gene_target         STRING  COMMENT 'CRISPR knocked-out gene (or non-targeting)',
        cell_type           STRING  COMMENT 'Predicted cell type',
        region              STRING  COMMENT 'Brain region',
        umap_1              DOUBLE  COMMENT 'UMAP embedding dim 1',
        umap_2              DOUBLE  COMMENT 'UMAP embedding dim 2',
        n_counts            INT     COMMENT 'Total UMI counts',
        batch               STRING  COMMENT 'Batch/sample identifier'
    """,
    "crispr_atlas_diff_expr": """
        gene_target         STRING  COMMENT 'CRISPR knocked-out gene',
        names               STRING  COMMENT 'Differentially expressed gene',
        group_name          STRING  COMMENT 'Cell type / cluster group',
        logfoldchanges      DOUBLE  COMMENT 'Log2 fold change (pert vs ctrl)',
        pvals_adj           DOUBLE  COMMENT 'BH-adjusted p-value',
        n_pert_matched      INT     COMMENT 'Matched perturbed cell count',
        n_ctrl_matched      INT     COMMENT 'Matched control cell count'
    """,
    "crispr_atlas_cell_metadata": """
        cell_barcode        STRING  COMMENT 'Unique cell identifier',
        gene_target         STRING  COMMENT 'CRISPR knocked-out gene',
        predicted_class     STRING  COMMENT 'Broad cell class (e.g. Glut, GABA)',
        predicted_subclass  STRING  COMMENT 'Subclass within class',
        neuron_type         STRING  COMMENT 'Neuron type label',
        region_level1       STRING  COMMENT 'Brain region (level 1)',
        region_level2       STRING  COMMENT 'Brain region (level 2)',
        passes_qc           BOOLEAN COMMENT 'Passes QC filter'
    """,
    "crispr_atlas_gene_metadata": """
        gene_target         STRING  COMMENT 'Knocked-out gene symbol',
        gene_id             STRING  COMMENT 'Ensembl gene ID',
        chromosome          STRING  COMMENT 'Chromosome',
        n_guides            INT     COMMENT 'Number of sgRNA guides',
        n_cells_targeted    INT     COMMENT 'Total cells with this KO'
    """,
    "crispr_atlas_ndd_genes": """
        gene                STRING  COMMENT 'Gene symbol',
        gene_id             STRING  COMMENT 'Ensembl gene ID',
        chromosome          STRING  COMMENT 'Chromosome',
        FDR_TADA_ASD        DOUBLE  COMMENT 'FDR from TADA model (ASD)',
        FDR_TADA_DD         DOUBLE  COMMENT 'FDR from TADA model (DD)',
        FDR_TADA_NDD        DOUBLE  COMMENT 'FDR from TADA model (NDD)',
        p_TADA_ASD          DOUBLE  COMMENT 'TADA p-value (ASD)',
        p_TADA_DD           DOUBLE  COMMENT 'TADA p-value (DD)',
        p_TADA_NDD          DOUBLE  COMMENT 'TADA p-value (NDD)',
        ASD72               DOUBLE  COMMENT 'ASD-72 gene list score',
        DD309               DOUBLE  COMMENT 'DD-309 gene list score',
        NDD373              DOUBLE  COMMENT 'NDD-373 gene list score',
        SCZ244              DOUBLE  COMMENT 'SCZ-244 gene list score'
    """,
    "gene_expression_matrix": """
        gene_symbol         STRING  COMMENT 'Gene symbol (UPPER CASE)',
        cohort              STRING  COMMENT 'Study cohort (BioFINDER, ROSMAP)',
        tissue              STRING  COMMENT 'Tissue or brain region',
        sample_id           STRING  COMMENT 'Sample identifier',
        tpm                 DOUBLE  COMMENT 'Transcripts per million'
    """,
    "tahoe_100m": """
        cell_id             STRING  COMMENT 'Unique cell identifier',
        drug                STRING  COMMENT 'Drug compound name',
        moa                 STRING  COMMENT 'Mechanism of action',
        cell_line_id        STRING  COMMENT 'Cell line identifier',
        dose_um             DOUBLE  COMMENT 'Drug concentration (µM)',
        timepoint_h         INT     COMMENT 'Treatment timepoint (hours)',
        n_counts            INT     COMMENT 'Total UMI counts'
    """,
    "tahoe_100m_clustered": """
        cell_id             STRING  COMMENT 'Unique cell identifier',
        drug                STRING  COMMENT 'Drug compound name',
        moa                 STRING  COMMENT 'Mechanism of action',
        cell_line_id        STRING  COMMENT 'Cell line identifier',
        cluster             STRING  COMMENT 'Transcriptomic cluster assignment',
        umap_1              DOUBLE  COMMENT 'UMAP embedding dim 1',
        umap_2              DOUBLE  COMMENT 'UMAP embedding dim 2'
    """,
    "tahoe_100m_gene_vocab": """
        gene_index          INT     COMMENT 'Integer gene index in feature matrix',
        gene_id             STRING  COMMENT 'Ensembl gene ID',
        gene_name           STRING  COMMENT 'Gene symbol'
    """,
}

for tbl, cols in placeholder_ddls.items():
    fqn = f"{CFG.catalog}.{CFG.schema}.{tbl}"
    spark.sql(f"""
        CREATE TABLE IF NOT EXISTS {fqn} (
            {cols}
        )
        USING DELTA
        COMMENT '{_PLACEHOLDER_COMMENT}'
    """)
    print(f"  ✅ {tbl}")
print("\nAll placeholder tables created.")

# COMMAND ----------

# DBTITLE 1,Step 4b: Placeholder tables — neuroplex_* biomedical sources
# Schemas derived from SELECT queries in neuro_mcp_server/neuroplex_query.py
# All tables use VARIANT payload for flexible JSON storage; same access pattern.
_NP_COMMENT = "PLACEHOLDER — correct schema, data ingestion pending from external source."

neuroplex_ddls = {
    "neuroplex_opentargets": """
        gene_symbol   STRING  COMMENT 'Gene symbol (UPPER CASE)',
        disease       STRING  COMMENT 'Disease name',
        payload       VARIANT COMMENT 'Full OpenTargets JSON (score, therapeuticAreas, etc.)'
    """,
    "neuroplex_gnomad": """
        gene_symbol   STRING  COMMENT 'Gene symbol (UPPER CASE)',
        title         STRING  COMMENT 'Record title (variant HGVS or Constraint)',
        payload       VARIANT COMMENT 'Full gnomAD JSON (clinical_significance, constraint, allele freq, etc.)'
    """,
    "neuroplex_kegg": """
        gene_symbol   STRING  COMMENT 'Gene symbol (UPPER CASE)',
        title         STRING  COMMENT 'Pathway entry title',
        summary       STRING  COMMENT 'Pathway description',
        payload       VARIANT COMMENT 'Full KEGG JSON (pathway_id, pathway_name, etc.)'
    """,
    "neuroplex_cbioportal": """
        gene_symbol   STRING  COMMENT 'Gene symbol (UPPER CASE)',
        title         STRING  COMMENT 'Mutation entry title',
        disease       STRING  COMMENT 'Cancer type / disease',
        summary       STRING  COMMENT 'Mutation summary',
        payload       VARIANT COMMENT 'Full cBioPortal JSON (mutation_type, protein_change, cancer_study, etc.)'
    """,
    "neuroplex_uniprot": """
        gene_symbol   STRING  COMMENT 'Gene symbol (UPPER CASE)',
        title         STRING  COMMENT 'UniProt entry title',
        summary       STRING  COMMENT 'Protein function summary',
        payload       VARIANT COMMENT 'Full UniProt JSON (protein_name, function, subcellular_location, etc.)'
    """,
    "neuroplex_ncbi_gtr": """
        gene_symbol   STRING  COMMENT 'Gene symbol (UPPER CASE)',
        title         STRING  COMMENT 'Test name',
        summary       STRING  COMMENT 'Test summary',
        disease       STRING  COMMENT 'Disease indication',
        payload       VARIANT COMMENT 'Full NCBI GTR JSON (test_type, lab_name, etc.)'
    """,
}

for tbl, cols in neuroplex_ddls.items():
    fqn = f"{CFG.catalog}.{CFG.schema}.{tbl}"
    spark.sql(f"""
        CREATE TABLE IF NOT EXISTS {fqn} (
            {cols}
        )
        USING DELTA
        COMMENT '{_NP_COMMENT}'
    """)
    print(f"  \u2705 {tbl}")
print("\nAll neuroplex_* placeholder tables created.")

# COMMAND ----------

# DBTITLE 1,Step 5: Grant SELECT to app service principal
SP_APP_ID = "79e31bca-0a8a-4fa8-b858-976114fb1f9b"  # neuro-crispr-mcp SP
tables = [
    "crisprbrain_screens",
    "gtex_brain_expression",
    "lincs_l1000_signatures",
    "chembl_orexin_pharmacology",
    "neuroplex_query_log",
    "wholebrain_crispr_atlas",
    "crispr_atlas_diff_expr",
    "crispr_atlas_cell_metadata",
    "crispr_atlas_gene_metadata",
    "crispr_atlas_ndd_genes",
    "gene_expression_matrix",
    "tahoe_100m",
    "tahoe_100m_clustered",
    "tahoe_100m_gene_vocab",
    # neuroplex_* biomedical sources
    "neuroplex_opentargets",
    "neuroplex_gnomad",
    "neuroplex_kegg",
    "neuroplex_cbioportal",
    "neuroplex_uniprot",
    "neuroplex_ncbi_gtr",
]
for t in tables:
    try:
        spark.sql(f"GRANT SELECT ON TABLE {CFG.catalog}.{CFG.schema}.{t} TO `{SP_APP_ID}`")
        print(f"  ✅ {t}")
    except Exception as e:
        print(f"  ⚠️  {t}: {e}")

# COMMAND ----------

# DBTITLE 1,Final status
expected = set(tables) | {"neuroplex_query_log"}
actual = {r.tableName for r in spark.sql(f"SHOW TABLES IN {CFG.catalog}.{CFG.schema}").collect()}
print(f"{'Table':<45} {'Status'}")
print("-" * 55)
for t in sorted(expected):
    status = "✅ EXISTS" if t in actual else "❌ MISSING"
    print(f"  {t:<43} {status}")
print(f"\n{len(expected & actual)}/{len(expected)} tables present")
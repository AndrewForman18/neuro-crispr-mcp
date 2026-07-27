# Databricks notebook source
# /// script
# [tool.databricks.environment]
# environment_version = "5"
# ///
# DBTITLE 1,Environment Setup
# MAGIC %md
# MAGIC # NeuroPlex — Deploy All Prod Tables
# MAGIC Orchestrates ingestion of all NeuroPlex tables into `dhbl_discovery_us_prod.genesis_schema`.
# MAGIC Runs: GTEx / LINCS / ChEMBL → CRISPRbrain screens → 22 neuroplex_* registry tables → GRANT permissions.

# COMMAND ----------

# DBTITLE 1,Step 1 — Pharmacology & Expression (GTEx / LINCS L1000 / ChEMBL)
import os, sys
from pathlib import PurePosixPath

# Point all ingestion code at prod
os.environ['NEUROPLEX_ENV']     = 'prod'
os.environ['NEUROPLEX_CATALOG'] = 'dhbl_discovery_us_prod'
os.environ['NEUROPLEX_SCHEMA']  = 'genesis_schema'

# Resolve repo root from this notebook's path
ctx = dbutils.notebook.entry_point.getDbutils().notebook().getContext()
repo_root = str(PurePosixPath(ctx.notebookPath().get()).parent)
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

from config.neuroplex_config import load_config
cfg = load_config()
print(f"✅ env={cfg.environment}  catalog={cfg.catalog}  schema={cfg.schema}")
print(f"   warehouse={cfg.sql_warehouse_id}")

# COMMAND ----------

# DBTITLE 1,Step 2 — CRISPRbrain Screens
# MAGIC %md
# MAGIC ## Step 1 — Pharmacology & Expression (GTEx, LINCS L1000, ChEMBL)

# COMMAND ----------

# DBTITLE 1,Step 3 — 22 NeuroPlex Registry Tables
# MAGIC %run ./02_ingest_pharmacology

# COMMAND ----------

# DBTITLE 1,Grant SELECT to neuro-crispr-mcp Service Principal
# MAGIC %md
# MAGIC ## Step 2 — CRISPRbrain Screens

# COMMAND ----------

# DBTITLE 1,Table Inventory
# MAGIC %run ./01_ingest_crisprbrain

# COMMAND ----------

# DBTITLE 1,Step 3 — 22 NeuroPlex Registry Tables
# MAGIC %md
# MAGIC ## Step 3 — 22 NeuroPlex Registry Tables (OpenTargets, gnomAD, ClinVar, etc.)

# COMMAND ----------

# DBTITLE 1,Run NeuroPlex Registry Ingestion
# MAGIC %run ./ingestion/03_run_ingestion_pipeline

# COMMAND ----------

# DBTITLE 1,Step 4 — Grant SELECT
# MAGIC %md
# MAGIC ## Step 4 — Grant SELECT to neuro-crispr-mcp Service Principal

# COMMAND ----------

# DBTITLE 1,Grant SELECT to neuro-crispr-mcp Service Principal
PROD_SP = "2a8be974-109e-48a9-8d86-06e567ac7c3a"
CATALOG  = "dhbl_discovery_us_prod"
SCHEMA   = "genesis_schema"

tables = spark.sql(f"SHOW TABLES IN {CATALOG}.{SCHEMA}").collect()
granted, failed = [], []
for row in tables:
    tbl = row["tableName"]
    try:
        spark.sql(f"GRANT SELECT ON TABLE {CATALOG}.{SCHEMA}.{tbl} TO `{PROD_SP}`")
        granted.append(tbl)
    except Exception as e:
        failed.append((tbl, str(e)))

print(f"✅ Granted SELECT on {len(granted)} tables to SP {PROD_SP}")
if failed:
    print(f"⚠️  Failed on {len(failed)} tables:")
    for t, e in failed: print(f"   {t}: {e}")

# COMMAND ----------

# DBTITLE 1,Step 5 — Verify Table Coverage
# MAGIC %md
# MAGIC ## Step 5 — Verify Table Coverage

# COMMAND ----------

# DBTITLE 1,Table Inventory
# MAGIC %sql
# MAGIC SELECT
# MAGIC   table_name,
# MAGIC   CAST(COALESCE(row_count, -1) AS BIGINT) AS row_count,
# MAGIC   CAST(COALESCE(size_bytes / 1024 / 1024, 0) AS BIGINT) AS size_mb
# MAGIC FROM dhbl_discovery_us_prod.information_schema.tables
# MAGIC WHERE table_schema = 'genesis_schema'
# MAGIC ORDER BY table_name

# COMMAND ----------

# DBTITLE 1,Step 3 — 22 NeuroPlex Registry Tables
# MAGIC %md
# MAGIC ## Step 3 — 22 NeuroPlex Registry Tables (OpenTargets, gnomAD, ClinVar, etc.)

# COMMAND ----------

# DBTITLE 1,Run NeuroPlex Registry Ingestion
# MAGIC %run ./ingestion/03_run_ingestion_pipeline

# COMMAND ----------

# DBTITLE 1,Run NeuroPlex Registry Ingestion
# MAGIC %run ./ingestion/03_run_ingestion_pipeline

# COMMAND ----------

# DBTITLE 1,Step 4 — Grant SELECT
# MAGIC %md
# MAGIC ## Step 4 — Grant SELECT to neuro-crispr-mcp Service Principal

# COMMAND ----------

# DBTITLE 1,Step 4 — Grant SELECT
# MAGIC %md
# MAGIC ## Step 4 — Grant SELECT to neuro-crispr-mcp Service Principal

# COMMAND ----------

# DBTITLE 1,Grant SELECT to neuro-crispr-mcp Service Principal
PROD_SP = "2a8be974-109e-48a9-8d86-06e567ac7c3a"
CATALOG  = "dhbl_discovery_us_prod"
SCHEMA   = "genesis_schema"

tables = spark.sql(f"SHOW TABLES IN {CATALOG}.{SCHEMA}").collect()
granted, failed = [], []
for row in tables:
    tbl = row["tableName"]
    try:
        spark.sql(f"GRANT SELECT ON TABLE {CATALOG}.{SCHEMA}.{tbl} TO `{PROD_SP}`")
        granted.append(tbl)
    except Exception as e:
        failed.append((tbl, str(e)))

print(f"✅ Granted SELECT on {len(granted)} tables to SP {PROD_SP}")
if failed:
    print(f"⚠️  Failed on {len(failed)} tables:")
    for t, e in failed: print(f"   {t}: {e}")

# COMMAND ----------

# DBTITLE 1,Step 5 — Verify Table Coverage
# MAGIC %md
# MAGIC ## Step 5 — Verify Table Coverage

# COMMAND ----------

# DBTITLE 1,Table Inventory
# MAGIC %sql
# MAGIC SELECT
# MAGIC   table_name,
# MAGIC   CAST(COALESCE(row_count, -1) AS BIGINT) AS row_count,
# MAGIC   CAST(COALESCE(size_bytes / 1024 / 1024, 0) AS BIGINT) AS size_mb
# MAGIC FROM dhbl_discovery_us_prod.information_schema.tables
# MAGIC WHERE table_schema = 'genesis_schema'
# MAGIC ORDER BY table_name

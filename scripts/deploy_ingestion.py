"""NeuroPlex Deployment Ingestion Script.

Run this script in a Databricks notebook cell after cloning the repo to a new workspace.
It creates all required tables and populates them from public APIs.

Usage (in a notebook cell):
    %run ./scripts/deploy_ingestion

Or execute directly:
    exec(open('/Workspace/Users/<you>/neuro-crispr-mcp-repo/scripts/deploy_ingestion.py').read())

Prerequisites:
    - Unity Catalog with target catalog/schema created
    - config/neuroplex_env.yml populated for the target environment
    - NEUROPLEX_ENV set to the target environment name
    - pip install requests pyyaml
"""

import os
import sys
import json
import time
from pathlib import Path

# Resolve repo root (this script lives in scripts/)
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from config.neuroplex_config import load_config

CFG = load_config()
CATALOG = CFG.catalog
SCHEMA = CFG.schema

print(f"="*60)
print(f"NeuroPlex Deployment Ingestion")
print(f"="*60)
print(f"Environment: {CFG.environment}")
print(f"Target:      {CATALOG}.{SCHEMA}")
print(f"Warehouse:   {CFG.sql_warehouse_id}")
print(f"Repo root:   {REPO_ROOT}")
print()

# ═══════════════════════════════════════════════════════════════════════════════
# Phase 1: Create Schema & Tables
# ═══════════════════════════════════════════════════════════════════════════════

print("[Phase 1] Creating schema and tables...")

CREATE_SCHEMA_SQL = f"CREATE SCHEMA IF NOT EXISTS {CATALOG}.{SCHEMA}"

# Common table schema for neuroplex_* tables
NEUROPLEX_TABLE_DDL = """
CREATE TABLE IF NOT EXISTS {catalog}.{schema}.{table} (
    record_id STRING NOT NULL,
    source_key STRING NOT NULL,
    gene_symbol STRING,
    disease STRING,
    drug STRING,
    title STRING,
    summary STRING,
    payload STRING,
    ingested_at TIMESTAMP,
    source_updated_at STRING
)
USING DELTA
COMMENT 'NeuroPlex {source_name} data (auto-created by deploy_ingestion.py)'
"""

# Core app tables (non-neuroplex)
CORE_TABLES_DDL = [
    f"""
    CREATE TABLE IF NOT EXISTS {CATALOG}.{SCHEMA}.crisprbrain_screens (
        gene STRING,
        screen_name STRING,
        cell_type STRING,
        crispr_mode STRING,
        phenotype_name STRING,
        genotype STRING,
        log2fc DOUBLE,
        pvalue DOUBLE,
        fdr DOUBLE,
        phenotype_score DOUBLE,
        source STRING,
        source_file STRING,
        sheet_name STRING
    ) USING DELTA
    COMMENT 'CRISPRbrain human iPSC CRISPR screens'
    """,
    f"""
    CREATE TABLE IF NOT EXISTS {CATALOG}.{SCHEMA}.neuroplex_query_log (
        query_id STRING,
        user_email STRING,
        question STRING,
        datasets STRING,
        response_summary STRING,
        tools_called STRING,
        latency_ms BIGINT,
        created_at TIMESTAMP
    ) USING DELTA
    COMMENT 'NeuroPlex query audit log'
    """,
]

# NeuroPlex source tables
NEUROPLEX_SOURCES = [
    ("neuroplex_opentargets", "OpenTargets disease associations"),
    ("neuroplex_gnomad", "gnomAD v4 variants and gene constraint"),
    ("neuroplex_kegg", "KEGG pathway memberships"),
    ("neuroplex_cbioportal", "cBioPortal somatic mutations"),
    ("neuroplex_uniprot", "UniProt protein annotations"),
    ("neuroplex_ncbi_gtr", "NCBI Genetic Testing Registry"),
    ("neuroplex_civic", "CIViC clinical variant evidence"),
    ("neuroplex_reactome", "Reactome pathway data"),
    ("neuroplex_monarch", "Monarch Initiative phenotype data"),
    ("neuroplex_hpa", "Human Protein Atlas expression"),
]


def run_ddl(spark_session):
    """Create all tables. Requires an active SparkSession."""
    spark_session.sql(CREATE_SCHEMA_SQL)
    print(f"  Schema {CATALOG}.{SCHEMA} ensured.")

    for table_name, source_name in NEUROPLEX_SOURCES:
        ddl = NEUROPLEX_TABLE_DDL.format(
            catalog=CATALOG, schema=SCHEMA,
            table=table_name, source_name=source_name
        )
        spark_session.sql(ddl)
        print(f"  Table {table_name} ensured.")

    for ddl in CORE_TABLES_DDL:
        spark_session.sql(ddl)
    print(f"  Core tables ensured (crisprbrain_screens, neuroplex_query_log).")
    print()


# ═══════════════════════════════════════════════════════════════════════════════
# Phase 2: Ingest NeuroPlex Sources (OpenTargets, gnomAD, KEGG, etc.)
# ═══════════════════════════════════════════════════════════════════════════════

# Key neuroscience gene panel for ingestion
GENE_PANEL = [
    # Orexin system
    "HCRT", "HCRTR1", "HCRTR2",
    # Alzheimer's / FTD
    "PSEN1", "PSEN2", "APP", "MAPT", "GRN", "TARDBP", "C9orf72",
    # ALS
    "SOD1", "FUS", "OPTN", "TBK1",
    # NDD risk genes (top TADA hits)
    "CHD8", "SCN2A", "SYNGAP1", "DYRK1A", "ADNP",
    # Sleep/circadian
    "CLOCK", "PER2", "CRY1",
    # Parkinson's
    "LRRK2", "SNCA", "PARK7", "PINK1",
]


def run_neuroplex_ingestion(spark_session):
    """Run all NeuroPlex ingestors for the gene panel."""
    from ingestion.source_registry import SOURCE_MAP
    from ingestion.ingestors.opentargets import OpenTargetsIngestor
    from ingestion.ingestors.gnomad import GnomadIngestor
    from ingestion.ingestors.kegg import KeggIngestor
    from ingestion.ingestors.cbioportal import CbioportalIngestor
    from ingestion.ingestors.uniprot import UniprotIngestor
    from ingestion.ingestors.ncbi_gtr import NcbiGtrIngestor

    ingestors = [
        ("opentargets", OpenTargetsIngestor),
        ("gnomad", GnomadIngestor),
        ("kegg", KeggIngestor),
        ("cbioportal", CbioportalIngestor),
        ("uniprot", UniprotIngestor),
        ("ncbi_gtr", NcbiGtrIngestor),
    ]

    results = []
    for source_key, IngestorClass in ingestors:
        print(f"\n[{source_key}] Starting ingestion...")
        try:
            ingestor = IngestorClass(SOURCE_MAP[source_key])
            for gene in GENE_PANEL:
                try:
                    summary = ingestor.run(spark=spark_session, gene=gene, limit=50)
                    results.append(summary)
                    print(f"  {gene}: {summary.get('records_written', 0)} records")
                except Exception as e:
                    print(f"  {gene}: ERROR - {e}")
                    results.append({"source": source_key, "gene": gene, "error": str(e)})
        except Exception as e:
            print(f"  FAILED to initialize {source_key}: {e}")

    successful = [r for r in results if "error" not in r]
    failed = [r for r in results if "error" in r]
    print(f"\n[NeuroPlex] Complete: {len(successful)} successful, {len(failed)} failed")
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# Phase 3: Ingest CRISPRbrain Screens
# ═══════════════════════════════════════════════════════════════════════════════

def run_crisprbrain_ingestion(spark_session):
    """Ingest CRISPRbrain data from Cell 2022 supplementary files."""
    import requests
    import pandas as pd
    from io import BytesIO

    SUPP_URLS = [
        "https://ars.els-cdn.com/content/image/1-s2.0-S0092867422005979-mmc1.xlsx",
        "https://ars.els-cdn.com/content/image/1-s2.0-S0092867422005979-mmc2.xlsx",
        "https://ars.els-cdn.com/content/image/1-s2.0-S0092867422005979-mmc3.xlsx",
    ]

    all_dfs = []
    for url in SUPP_URLS:
        fname = url.split("/")[-1]
        print(f"  Downloading {fname}...")
        try:
            r = requests.get(url, timeout=120)
            r.raise_for_status()
            xls = pd.ExcelFile(BytesIO(r.content))
            for sheet in xls.sheet_names:
                try:
                    df = pd.read_excel(xls, sheet_name=sheet)
                    if "Gene" in df.columns or "gene" in df.columns:
                        df.columns = [c.strip() for c in df.columns]
                        if "gene" in df.columns:
                            df = df.rename(columns={"gene": "Gene"})
                        df["source_file"] = fname
                        df["sheet_name"] = sheet
                        all_dfs.append(df)
                except Exception:
                    pass
        except Exception as e:
            print(f"    ERROR downloading {fname}: {e}")

    if all_dfs:
        combined = pd.concat(all_dfs, ignore_index=True)
        print(f"  Total CRISPRbrain rows: {len(combined):,}")

        # Write to Delta
        target = f"{CATALOG}.{SCHEMA}.crisprbrain_screens"
        sdf = spark_session.createDataFrame(combined.astype(str))
        sdf.write.format("delta").mode("overwrite").saveAsTable(target)
        print(f"  Written to {target}")
    else:
        print("  WARNING: No CRISPRbrain data downloaded.")


# ═══════════════════════════════════════════════════════════════════════════════
# Main Execution
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    """Full deployment ingestion pipeline."""
    from pyspark.sql import SparkSession
    spark = SparkSession.getActiveSession()
    if not spark:
        print("ERROR: No active SparkSession. Run this in a Databricks notebook.")
        return

    start = time.time()

    # Phase 1: DDL
    print("[Phase 1] Creating tables...")
    run_ddl(spark)

    # Phase 2: NeuroPlex sources
    print("\n[Phase 2] Ingesting NeuroPlex sources...")
    run_neuroplex_ingestion(spark)

    # Phase 3: CRISPRbrain
    print("\n[Phase 3] Ingesting CRISPRbrain...")
    run_crisprbrain_ingestion(spark)

    elapsed = time.time() - start
    print(f"\n{'='*60}")
    print(f"Deployment ingestion complete in {elapsed/60:.1f} minutes.")
    print(f"{'='*60}")
    print()
    print("NOTE: The following large datasets require separate ingestion:")
    print("  - wholebrain_crispr_atlas (7.7M cells) — contact platform team for Parquet source")
    print("  - crispr_atlas_diff_expr (745M rows) — contact platform team for Parquet source")
    print("  - gene_expression_matrix (21M rows) — BioFINDER/ROSMAP cohort data")
    print("  - gtex_brain_expression — GTEx v8 brain tissue download")
    print("  - lincs_l1000_signatures — LINCS L1000 bulk download")
    print("  - chembl_orexin_pharmacology — ChEMBL API extraction")
    print("  - tahoe_100m (95M rows) — Tahoe-100M drug atlas Parquet")
    print()
    print("Run 01_ingest_crisprbrain and 02_ingest_pharmacology notebooks")
    print("for the remaining datasets after this script completes.")


if __name__ == "__main__":
    main()

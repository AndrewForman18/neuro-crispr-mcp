# Databricks notebook source
# /// script
# [tool.databricks.environment]
# environment_version = "5"
# ///
# DBTITLE 1,Install — huggingface_hub + pyarrow
# MAGIC %pip install -q huggingface_hub pyarrow --upgrade

# COMMAND ----------

# DBTITLE 1,Config — environment, tables, HF repo
import os, sys, json, time, logging
import pyarrow.parquet as pq
import pandas as pd
from pyspark.sql import functions as F
from pyspark.sql.types import StructType, StructField, StringType, DoubleType, IntegerType, BooleanType

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("crispr_atlas_etl")

os.environ["NEUROPLEX_ENV"] = "val"
sys.path.insert(0, "/Workspace/Users/andrew_forman@eisai.com/neuro-crispr-mcp")
for mod in [k for k in sys.modules if k.startswith("config")]:
    del sys.modules[mod]
from config.neuroplex_config import load_config

CFG = load_config()
CAT, SCH = CFG.catalog, CFG.schema

# Target Delta tables
TABLE_ATLAS     = f"{CAT}.{SCH}.wholebrain_crispr_atlas"
TABLE_DIFF_EXPR = f"{CAT}.{SCH}.crispr_atlas_diff_expr"
TABLE_CELL_META = f"{CAT}.{SCH}.crispr_atlas_cell_metadata"
TABLE_GENE_META = f"{CAT}.{SCH}.crispr_atlas_gene_metadata"
TABLE_NDD_GENES = f"{CAT}.{SCH}.crispr_atlas_ndd_genes"

# HuggingFace source
CRISPR_REPO  = "perturbai/wholebrain_crispr_atlas"
HF_REPO_TYPE = "dataset"
HF_TOKEN     = os.environ.get("HUGGING_FACE_HUB_TOKEN", None)

# HF cache goes to /tmp — lost on cluster restart but avoids re-downloading within a session
# For persistent caching, set to a Volume path:
# os.environ["HF_HOME"] = "/Volumes/dhbl_discovery_us_val/genesis_schema/data/hf_cache"

print(f"\u2705 Environment : {CFG.environment}")
for name, tbl in [("Atlas",TABLE_ATLAS),("DiffExpr",TABLE_DIFF_EXPR),("CellMeta",TABLE_CELL_META),
                   ("GeneMeta",TABLE_GENE_META),("NddGenes",TABLE_NDD_GENES)]:
    print(f"   {name:<12s}: {tbl}")
print(f"   HF repo   : {CRISPR_REPO}  (token={'set' if HF_TOKEN else 'not set (anon OK)'})")

# COMMAND ----------

# DBTITLE 1,Discover — list HF files, print schema of cell data + DE parquet
# Reads ONLY schema metadata (not data rows) to avoid OOM.
# The DE parquet has ~745M rows — never load into pandas; use Spark only.

from huggingface_hub import list_repo_files, hf_hub_download

files_crispr = list(list_repo_files(CRISPR_REPO, repo_type=HF_REPO_TYPE, token=HF_TOKEN))

# Categorise
data_shards = sorted([f for f in files_crispr if f.startswith("data/") and f.endswith(".parquet")])
de_file     = next(f for f in files_crispr if "wilcoxon_de_results" in f)
ndd_csvs    = sorted([f for f in files_crispr if "ndd_perturbation" in f
                      and f.split(".")[-1] in ("csv", "tsv")])
h5ad_files  = sorted([f for f in files_crispr if f.endswith(".h5ad")])

print(f"Cell data shards : {len(data_shards)}")
print(f"DE results file  : {de_file}")
print(f"NDD CSV/TSV files : {len(ndd_csvs)}")
print(f"h5ad files        : {len(h5ad_files)}")

# ---- Cell data shard: schema only ----
print(f"\n--- Cell data schema (schema only, no row load) ---")
p_cell = hf_hub_download(CRISPR_REPO, data_shards[0], repo_type=HF_REPO_TYPE, token=HF_TOKEN)
pf_cell = pq.ParquetFile(p_cell)
print(f"  Rows per shard : {pf_cell.metadata.num_rows:,}  |  Row groups: {pf_cell.metadata.num_row_groups}")
print(f"  Fields ({len(pf_cell.schema_arrow)}):")
for field in pf_cell.schema_arrow:
    print(f"    {field.name:<40s}  {str(field.type)}")
# Show 3 rows with safe columns (no large embedding fields)
safe_cols = [f.name for f in pf_cell.schema_arrow if str(f.type) not in ("large_list<item: double>", "list<item: double>")][:12]
samp_cell = pf_cell.read_row_group(0, columns=safe_cols).to_pandas().head(3)
print(samp_cell.to_string())
del pf_cell

# ---- DE results: schema only (no row load) ----
print(f"\n--- DE results schema (schema only) ---")
p_de = hf_hub_download(CRISPR_REPO, de_file, repo_type=HF_REPO_TYPE, token=HF_TOKEN)
pf_de = pq.ParquetFile(p_de)
print(f"  Total rows : {pf_de.metadata.num_rows:,}  |  Row groups: {pf_de.metadata.num_row_groups}")
print(f"  Fields: {list(pf_de.schema_arrow.names)}")
# First 5 rows from row group 0
de_samp = pf_de.read_row_group(0).to_pandas().head(5)
print(de_samp.to_string())
del pf_de

print("\n--- NDD CSV file list ---")
for f in ndd_csvs:
    print(f"  {f}")

# COMMAND ----------

# DBTITLE 1,Build crispr_atlas_diff_expr — Spark reads DE parquet (745M rows)
# Confirmed schema (745,758,369 rows, 10.8 GB parquet):
# Names match server.py exactly — NO column renames needed.
#   names, scores, logfoldchanges, pvals, pvals_adj, group_name,
#   gene_target, control_label, n_pert_matched, n_ctrl_matched
#
# p_de was downloaded and cached in probe cell (cell 3).
# If kernel was restarted, re-run cell 3 first to re-download p_de.
#
# Expected runtime: 20-40 min (745M rows, shuffle to Delta).

# Serverless Spark Connect blocks file:// paths (/tmp/).
# Copy the cached DE file to a Volume first, then Spark reads from the Volume.
import shutil, os

VOL_DIR    = "/Volumes/dhbl_discovery_us_val/genesis_schema/data"
VOL_DE     = f"{VOL_DIR}/crispr_de_results.parquet"

if not os.path.exists(VOL_DE):
    size_gb = os.path.getsize(p_de) / 1e9
    print(f"Copying DE parquet ({size_gb:.1f} GB) from /tmp/ → Volume ...")
    t_copy = time.time()
    shutil.copy2(p_de, VOL_DE)
    print(f"  Copied in {(time.time()-t_copy)/60:.1f} min")
else:
    print(f"  DE parquet already in Volume: {VOL_DE}")

print(f"\nReading DE results with Spark from: {VOL_DE}")
print("  Expected: 745,758,369 rows — write will take ~20-40 min")

df_de_raw = spark.read.parquet(VOL_DE)
df_de_raw.show(5, truncate=80)

# Select the 7 columns server.py queries + keep scores/pvals/control_label as extras
DE_KEEP = ["gene_target", "names", "group_name", "logfoldchanges",
           "pvals_adj", "n_pert_matched", "n_ctrl_matched",
           "scores", "pvals", "control_label"]
final_de_cols = [c for c in DE_KEEP if c in df_de_raw.columns]
df_de_final = df_de_raw.select(*final_de_cols)

print(f"\n  Final columns : {final_de_cols}")
print(f"  Writing to {TABLE_DIFF_EXPR} ...")
t0 = time.time()
df_de_final.write.format("delta").mode("overwrite").option("overwriteSchema","true") \
    .partitionBy("gene_target") \
    .saveAsTable(TABLE_DIFF_EXPR)
print(f"  \u2705 Done in {(time.time()-t0)/60:.1f} min → {TABLE_DIFF_EXPR}")

# COMMAND ----------

# DBTITLE 1,Build wholebrain_crispr_atlas + crispr_atlas_cell_metadata
# Serverless Spark Connect blocks spark.read.parquet("file:///tmp/...").
# Workaround: use pyarrow to read each shard from /tmp/ (pure Python),
# then spark.createDataFrame(pandas_df) to lift into Spark — no file I/O via Spark.
# Per-shard pandas footprint: ~4-6 MB (9 string/bool cols × 25k rows, skipping expression arrays).
# Strategy: overwrite on first shard, append on remainder.

META_COLS_READ = [
    "cell_id", "gene_target", "predicted_class", "predicted_subclass",
    "neuron_type", "region_level1", "region_level2", "passes_qc",
    "batch", "num_rna_umi",
]

# Discover which requested columns actually exist in the parquet
p0 = hf_hub_download(CRISPR_REPO, data_shards[0], repo_type=HF_REPO_TYPE, token=HF_TOKEN)
available = set(pq.ParquetFile(p0).schema_arrow.names)
read_cols = [c for c in META_COLS_READ if c in available]
print(f"Reading columns per shard : {read_cols}")
print(f"Processing {len(data_shards)} shards ...\n")

atlas_mode = "overwrite"
meta_mode  = "overwrite"
t0 = time.time()

for i, shard in enumerate(data_shards):
    lp = hf_hub_download(CRISPR_REPO, shard, repo_type=HF_REPO_TYPE, token=HF_TOKEN)

    # pyarrow reads /tmp/ path — returns only the 10 metadata columns
    pa_table = pq.read_table(lp, columns=read_cols)
    pd_df    = pa_table.to_pandas()
    del pa_table

    df_shard = spark.createDataFrame(pd_df)   # driver → Spark, no file I/O
    del pd_df

    # ---- wholebrain_crispr_atlas ----
    df_atlas_s = df_shard.select(
        F.col("cell_id").alias("cell_barcode"),
        F.col("gene_target"),
        F.col("predicted_class").alias("cell_type"),
        F.col("region_level1").alias("region"),
        F.lit(None).cast("double").alias("umap_1"),
        F.lit(None).cast("double").alias("umap_2"),
        F.col("num_rna_umi").cast("integer").alias("n_counts"),
        F.col("batch"),
    )
    df_atlas_s.write.format("delta").mode(atlas_mode) \
        .option("overwriteSchema", str(atlas_mode == "overwrite").lower()) \
        .saveAsTable(TABLE_ATLAS)
    atlas_mode = "append"

    # ---- crispr_atlas_cell_metadata ----
    df_meta_s = df_shard.select(
        F.col("cell_id").alias("cell_barcode"),
        F.col("gene_target"),
        F.col("predicted_class"),
        F.col("predicted_subclass"),
        F.col("neuron_type"),
        F.col("region_level1"),
        F.col("region_level2"),
        F.col("passes_qc"),
        F.col("batch"),
    )
    df_meta_s.write.format("delta").mode(meta_mode) \
        .option("overwriteSchema", str(meta_mode == "overwrite").lower()) \
        .saveAsTable(TABLE_CELL_META)
    meta_mode = "append"

    if (i + 1) % 30 == 0 or i == len(data_shards) - 1:
        elapsed = time.time() - t0
        rate    = (i + 1) / elapsed
        eta     = (len(data_shards) - i - 1) / rate if rate > 0 else 0
        print(f"  {i+1:>3}/{len(data_shards)}  {elapsed:>6.0f}s elapsed  ~{eta:.0f}s remaining")

print(f"\n\u2705 All {len(data_shards)} shards written in {(time.time()-t0)/60:.1f} min")

# COMMAND ----------

# DBTITLE 1,Build crispr_atlas_gene_metadata — per-KO target summary
# Derive gene metadata from cell_metadata:
# n_guides comes from the DE file (unique sgRNAs), n_cells_targeted from cell counts.
# Ensembl IDs + chromosome from a BioMart lookup (optional) or left empty.

from pyspark.sql import Window

# Derive from cell metadata
df_gene_meta = (
    spark.table(TABLE_CELL_META)
    .filter(F.col("gene_target").isNotNull())
    .filter(F.col("gene_target") != "non-targeting")
    .groupBy("gene_target")
    .agg(
        F.count("*").alias("n_cells_targeted"),
    )
    .withColumn("gene_id",    F.lit(None).cast("string"))   # fill via BioMart if needed
    .withColumn("chromosome", F.lit(None).cast("string"))
    .withColumn("n_guides",   F.lit(None).cast("int"))
    .select("gene_target", "gene_id", "chromosome", "n_guides", "n_cells_targeted")
)

# Add n_guides from DE file only if a guide-identity column exists
de_cols = set(spark.table(TABLE_DIFF_EXPR).columns)
guide_col = next((c for c in ["sgRNA", "guide", "guide_id", "num_guides"] if c in de_cols), None)
if guide_col:
    df_guides = (
        spark.table(TABLE_DIFF_EXPR)
        .filter(F.col("gene_target").isNotNull())
        .groupBy("gene_target")
        .agg(F.countDistinct(guide_col).alias("n_guides"))
    )
    df_gene_meta = df_gene_meta.drop("n_guides").join(df_guides, on="gene_target", how="left")
else:
    print("  No guide-identity column in DE table — deriving n_guides from cell_metadata")
    df_guides2 = (
        spark.table(TABLE_CELL_META)
        .filter(F.col("gene_target").isNotNull())
        .filter(F.col("gene_target") != "non-targeting")
        .groupBy("gene_target")
        .agg(F.countDistinct(F.col("batch")).alias("n_guides"))  # batches as proxy
    )
    df_gene_meta = df_gene_meta.drop("n_guides").join(df_guides2, on="gene_target", how="left")

df_gene_meta.write.format("delta").mode("overwrite").option("overwriteSchema","true").saveAsTable(TABLE_GENE_META)
print(f"  \u2705 {TABLE_GENE_META}  ({df_gene_meta.count():,} KO targets)")
display(df_gene_meta.orderBy("n_cells_targeted", ascending=False).limit(10))

# COMMAND ----------

# DBTITLE 1,Build crispr_atlas_ndd_genes — from NDD perturbation CSVs
# ASD_rare.csv is the NDD gene annotation table (18,128 genes genome-wide).
# It has EXACT server.py column names: gene, gene_id, chromosome,
# FDR_TADA_ASD/DD/NDD, p_TADA_ASD/DD/NDD, ASD72, DD309, NDD373, SCZ244.
# (GWS_results.csv is a DE results subset — wrong file for this table.)

for csv_f in ndd_csvs:
    if "ASD_rare" in csv_f:
        lp = hf_hub_download(CRISPR_REPO, csv_f, repo_type=HF_REPO_TYPE, token=HF_TOKEN)
        df_asd = pd.read_csv(lp)
        print(f"  ASD_rare.csv : {len(df_asd):,} rows")
        print(f"  Columns      : {list(df_asd.columns)}")
        break

# Select and rename to server.py schema (columns are already correctly named)
FINAL_COLS = ["gene", "gene_id", "chromosome",
              "FDR_TADA_ASD", "FDR_TADA_DD", "FDR_TADA_NDD",
              "p_TADA_ASD",   "p_TADA_DD",   "p_TADA_NDD",
              "ASD72", "DD309", "NDD373", "SCZ244"]
present = [c for c in FINAL_COLS if c in df_asd.columns]
missing = [c for c in FINAL_COLS if c not in df_asd.columns]
if missing:
    print(f"  Missing cols (will be NULL): {missing}")
    for c in missing:
        df_asd[c] = None

df_ndd_final = df_asd[FINAL_COLS].drop_duplicates(subset=["gene"])
print(f"\n  Genes        : {len(df_ndd_final):,}")
print(f"  ASD72 hits   : {(df_ndd_final['ASD72'] == 1).sum()}")
print(f"  DD309 hits   : {(df_ndd_final['DD309'] == 1).sum()}")
print(f"  NDD373 hits  : {(df_ndd_final['NDD373'] == 1).sum()}")
print(df_ndd_final[df_ndd_final["ASD72"] == 1].head(5).to_string())

df_ndd_spark = spark.createDataFrame(df_ndd_final.astype(str))
df_ndd_spark.write.format("delta").mode("overwrite").option("overwriteSchema","true").saveAsTable(TABLE_NDD_GENES)
print(f"  \u2705 {TABLE_NDD_GENES}  ({df_ndd_spark.count():,} genes)")

# COMMAND ----------

# DBTITLE 1,Grant SELECT to app SP on all five CRISPR Atlas tables
APP_SP = "79e31bca-0a8a-4fa8-b858-976114fb1f9b"  # neuro-crispr-mcp app service principal
CRISPR_TABLES = [TABLE_ATLAS, TABLE_DIFF_EXPR, TABLE_CELL_META, TABLE_GENE_META, TABLE_NDD_GENES]

for tbl in CRISPR_TABLES:
    try:
        spark.sql(f"GRANT SELECT ON TABLE {tbl} TO `{APP_SP}`")
        print(f"  \u2705 GRANT SELECT on {tbl.split('.')[-1]}")
    except Exception as e:
        print(f"  \u26a0\ufe0f  {tbl.split('.')[-1]}: {e}")

# COMMAND ----------

# DBTITLE 1,Validation — row counts across all five CRISPR Atlas tables
# MAGIC %sql
# MAGIC SELECT 'wholebrain_crispr_atlas'      AS tbl, COUNT(*) AS rows FROM dhbl_discovery_us_val.genesis_schema.wholebrain_crispr_atlas
# MAGIC UNION ALL
# MAGIC SELECT 'crispr_atlas_diff_expr'       AS tbl, COUNT(*) AS rows FROM dhbl_discovery_us_val.genesis_schema.crispr_atlas_diff_expr
# MAGIC UNION ALL
# MAGIC SELECT 'crispr_atlas_cell_metadata'   AS tbl, COUNT(*) AS rows FROM dhbl_discovery_us_val.genesis_schema.crispr_atlas_cell_metadata
# MAGIC UNION ALL
# MAGIC SELECT 'crispr_atlas_gene_metadata'   AS tbl, COUNT(*) AS rows FROM dhbl_discovery_us_val.genesis_schema.crispr_atlas_gene_metadata
# MAGIC UNION ALL
# MAGIC SELECT 'crispr_atlas_ndd_genes'       AS tbl, COUNT(*) AS rows FROM dhbl_discovery_us_val.genesis_schema.crispr_atlas_ndd_genes
# MAGIC ORDER BY tbl

# COMMAND ----------

# DBTITLE 1,Validation — DE results spot-check (Psen1 KO top effects)
# MAGIC %sql
# MAGIC -- Spot-check: top DE genes when Psen1 is knocked out
# MAGIC -- Expect significant effects in excitatory neurons (Glut, CTX cell types)
# MAGIC SELECT names AS affected_gene, group_name AS cell_type,
# MAGIC        ROUND(logfoldchanges, 3) AS log2fc, ROUND(pvals_adj, 6) AS padj,
# MAGIC        n_pert_matched, n_ctrl_matched
# MAGIC FROM dhbl_discovery_us_val.genesis_schema.crispr_atlas_diff_expr
# MAGIC WHERE gene_target = 'Psen1' AND pvals_adj < 0.05
# MAGIC ORDER BY ABS(logfoldchanges) DESC
# MAGIC LIMIT 20

# COMMAND ----------

# DBTITLE 1,Validation — NDD gene TADA scores sample
# MAGIC %sql
# MAGIC -- NDD genes with strongest ASD TADA evidence
# MAGIC SELECT gene, chromosome, FDR_TADA_ASD, FDR_TADA_DD, FDR_TADA_NDD, ASD72, DD309
# MAGIC FROM dhbl_discovery_us_val.genesis_schema.crispr_atlas_ndd_genes
# MAGIC WHERE FDR_TADA_ASD IS NOT NULL
# MAGIC ORDER BY CAST(FDR_TADA_ASD AS DOUBLE) ASC
# MAGIC LIMIT 20

# COMMAND ----------

# DBTITLE 1,GBA1 patch — append GBA1 to gnomAD, UniProt, KEGG, and OpenTargets tables
# Runs all three neuroplex GBA1 patches inline from this notebook.
# No need to navigate to 03/04/05 — all API calls + Delta appends happen here.

import os, sys, json, time, requests
from pyspark.sql import functions as F
from pyspark.sql.types import StructType, StructField, StringType

os.environ["NEUROPLEX_ENV"] = "val"
sys.path.insert(0, "/Workspace/Users/andrew_forman@eisai.com/neuro-crispr-mcp")
for mod in [k for k in sys.modules if k.startswith("config")]:
    del sys.modules[mod]
from config.neuroplex_config import load_config
CFG = load_config()
CAT, SCH = CFG.catalog, CFG.schema

REQUEST_DELAY = 1.5   # seconds between API calls

def _append_records(records, table, str_cols):
    """Append a list of {gene_symbol, title, payload} dicts to a Delta table."""
    if not records:
        print(f"  ⚠️  0 records — nothing to append to {table.split('.')[-1]}")
        return 0
    rows = [
        tuple(r.get(c, "") for c in str_cols)
        + (json.dumps(r.get("payload", {}), default=str),)
        for r in records
    ]
    schema = StructType(
        [StructField(c, StringType(), True) for c in str_cols]
        + [StructField("payload_str", StringType(), True)]
    )
    df = spark.createDataFrame(rows, schema).select(
        *[F.col(c) for c in str_cols],
        F.expr("parse_json(payload_str)").alias("payload"),
    )
    n = df.count()
    df.write.format("delta").mode("append").saveAsTable(table)
    print(f"  ✅ Appended {n} GBA1 rows → {table.split('.')[-1]}")
    return n


# ============================================================
# 1.  gnomAD — constraint + ClinVar  → neuroplex_gnomad
# ============================================================
GNOMAD_URL = "https://gnomad.broadinstitute.org/api"

def _gql(query):
    r = requests.post(GNOMAD_URL, json={"query": query}, timeout=30)
    r.raise_for_status()
    return r.json()

def fetch_gnomad_constraint(symbol):
    q = f"""{{
  gene(gene_symbol: "{symbol}", reference_genome: GRCh38) {{
    gene_id name
    gnomad_constraint {{ pLI oe_lof oe_lof_upper mis_z obs_lof exp_lof
      obs_mis exp_mis obs_syn exp_syn lof_z }}
  }}
}}"""
    try:
        d = _gql(q)
        gene = (d.get("data") or {}).get("gene")
        if not gene:
            return None
        return {"gene_symbol": symbol, "title": f"Constraint: {symbol}",
                "payload": {"gene_id": gene.get("gene_id"),
                            "constraint": gene.get("gnomad_constraint") or {}}}
    except Exception as e:
        print(f"  gnomAD constraint error: {e}")
        return None

def fetch_gnomad_clinvar(symbol, max_vars=500):
    q = f"""{{
  gene(gene_symbol: "{symbol}", reference_genome: GRCh38) {{
    clinvar_variants {{
      variant_id hgvsp hgvsc clinical_significance gold_stars
      major_consequence in_gnomad
    }}
  }}
}}"""
    try:
        d = _gql(q)
        gene = (d.get("data") or {}).get("gene")
        if not gene:
            return []
        variants = gene.get("clinvar_variants") or []
        return [{"gene_symbol": symbol,
                 "title": v.get("hgvsp") or v.get("hgvsc") or v.get("variant_id", ""),
                 "payload": v}
                for v in variants[:max_vars]]
    except Exception as e:
        print(f"  gnomAD ClinVar error: {e}")
        return []

print("[1/3] gnomAD: GBA1 constraint + ClinVar")
cr = fetch_gnomad_constraint("GBA1")
if cr:
    c = cr["payload"]["constraint"]
    print(f"  pLI={c.get('pLI')}  LOEUF={c.get('oe_lof_upper')}")
time.sleep(REQUEST_DELAY)
vars_ = fetch_gnomad_clinvar("GBA1")
print(f"  {len(vars_)} ClinVar variants")
gnomad_records = ([cr] if cr else []) + vars_
_append_records(gnomad_records, f"{CAT}.{SCH}.neuroplex_gnomad", ["gene_symbol", "title"])


# ============================================================
# 2.  UniProt + KEGG  → neuroplex_uniprot, neuroplex_kegg
# ============================================================
time.sleep(REQUEST_DELAY)

def fetch_uniprot_gba1():
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {"query": 'gene_exact:"GBA1" AND organism_id:9606 AND reviewed:true',
              "format": "json", "fields": "accession,id,protein_name,gene_names,organism_name,"
              "cc_function,cc_subcellular_location,sequence", "size": 1}
    try:
        r = requests.get(url, params=params, timeout=20)
        r.raise_for_status()
        results = r.json().get("results", [])
        if not results:
            # Fallback: try GBA (synonym)
            params["query"] = 'gene_exact:"GBA" AND organism_id:9606 AND reviewed:true AND accession:P04062'
            r = requests.get(url, params=params, timeout=20)
            r.raise_for_status()
            results = r.json().get("results", [])
        if not results:
            return None
        hit = results[0]
        acc  = hit.get("primaryAccession", "")
        pname = (hit.get("proteinDescription") or {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
        fn   = " ".join(c.get("texts", [{}])[0].get("value", "") for c in hit.get("comments", []) if c.get("commentType") == "FUNCTION")
        sloc = " ".join(c.get("subcellularLocations", [{}])[0].get("location", {}).get("value", "") for c in hit.get("comments", []) if c.get("commentType") == "SUBCELLULAR LOCATION")
        return {"gene_symbol": "GBA1", "title": f"{acc} {pname}", "summary": pname,
                "payload": {"uniprot_acc": acc, "protein_name": pname,
                            "function": fn[:500], "subcellular_location": sloc}}
    except Exception as e:
        print(f"  UniProt error: {e}")
        return None

print("\n[2/3] UniProt: GBA1")
u = fetch_uniprot_gba1()
if u:
    print(f"  {u['payload']['uniprot_acc']}  {u['payload']['protein_name']}")
_append_records([u] if u else [], f"{CAT}.{SCH}.neuroplex_uniprot", ["gene_symbol", "title", "summary"])
time.sleep(REQUEST_DELAY)

# KEGG: fetch all human pathways once, then filter for GBA1
def fetch_kegg_gba1():
    try:
        r = requests.get("http://rest.kegg.jp/list/pathway/hsa", timeout=20)
        r.raise_for_status()
        all_pathways = {}
        for line in r.text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) == 2:
                all_pathways[parts[0].replace("path:", "")] = parts[1]
        # Find pathways containing GBA1/GBA
        r2 = requests.get("http://rest.kegg.jp/link/pathway/hsa:2629", timeout=20)  # GBA NCBI gene 2629
        r2.raise_for_status()
        pathway_ids = [line.split("\t")[1].replace("path:", "") for line in r2.text.strip().split("\n") if "\t" in line]
        return [{"gene_symbol": "GBA1",
                 "title": pid,
                 "summary": all_pathways.get(pid, pid),
                 "payload": {"pathway_id": pid, "pathway_name": all_pathways.get(pid, pid)}}
                for pid in pathway_ids if pid in all_pathways]
    except Exception as e:
        print(f"  KEGG error: {e}")
        return []

print("KEGG: GBA1")
kr = fetch_kegg_gba1()
print(f"  {len(kr)} KEGG pathways")
if kr:
    for p in kr[:5]:
        print(f"    {p['payload']['pathway_id']}: {p['payload']['pathway_name']}")
_append_records(kr, f"{CAT}.{SCH}.neuroplex_kegg", ["gene_symbol", "title", "summary"])
time.sleep(REQUEST_DELAY)


# ============================================================
# 3.  OpenTargets  → neuroplex_opentargets
# ============================================================
OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"

SEARCH_Q = """
query($q:String!){
  search(queryString:$q, entityNames:["target"]) {
    hits { id entity }
  }
}"""

ASSOC_Q = """
query($id:String!){
  target(ensemblId:$id){
    approvedSymbol
    associatedDiseases(page:{index:0,size:100},enableIndirect:true){
      rows {
        disease { id name therapeuticAreas { id name } }
        score
      }
    }
  }
}"""

def _ot_gql(query, variables):
    r = requests.post(OT_URL, json={"query": query, "variables": variables}, timeout=30)
    r.raise_for_status()
    return r.json()

print("\n[3/3] OpenTargets: GBA1")
try:
    # Resolve symbol → Ensembl ID
    sr = _ot_gql(SEARCH_Q, {"q": "GBA1"})
    hits = (sr.get("data") or {}).get("search", {}).get("hits", [])
    target_hits = [h for h in hits if h.get("entity") == "target"]
    if not target_hits:
        # Fallback to GBA
        sr = _ot_gql(SEARCH_Q, {"q": "GBA"})
        hits = (sr.get("data") or {}).get("search", {}).get("hits", [])
        target_hits = [h for h in hits if h.get("entity") == "target"]
    ensembl_id = target_hits[0]["id"] if target_hits else None
    print(f"  Ensembl ID: {ensembl_id}")
except Exception as e:
    print(f"  OT search error: {e}")
    ensembl_id = None

ot_records = []
if ensembl_id:
    time.sleep(REQUEST_DELAY)
    try:
        ar = _ot_gql(ASSOC_Q, {"id": ensembl_id})
        target_data = (ar.get("data") or {}).get("target", {})
        rows = (target_data.get("associatedDiseases") or {}).get("rows", [])
        for row in rows:
            dis = row.get("disease", {})
            areas = dis.get("therapeuticAreas", [{}])
            ot_records.append({
                "gene_symbol": "GBA1",
                "disease": dis.get("name", ""),
                "payload": {"score": row.get("score"),
                            "disease_id": dis.get("id"),
                            "therapeuticAreas": [{"id": a.get("id"), "name": a.get("name")} for a in areas]}
            })
        print(f"  {len(ot_records)} disease associations")
        if ot_records:
            top = max(ot_records, key=lambda r: r["payload"].get("score") or 0)
            print(f"  Top: {top['disease']}  score={top['payload']['score']:.3f}")
    except Exception as e:
        print(f"  OT associations error: {e}")

# OT table uses different str_cols: gene_symbol, disease
if ot_records:
    rows_ot = [
        (r["gene_symbol"], r["disease"], json.dumps(r["payload"], default=str))
        for r in ot_records
    ]
    schema_ot = StructType([
        StructField("gene_symbol", StringType(), True),
        StructField("disease",     StringType(), True),
        StructField("payload_str", StringType(), True),
    ])
    df_ot = spark.createDataFrame(rows_ot, schema_ot).select(
        F.col("gene_symbol"), F.col("disease"),
        F.expr("parse_json(payload_str)").alias("payload"),
    )
    n_ot = df_ot.count()
    df_ot.write.format("delta").mode("append").saveAsTable(f"{CAT}.{SCH}.neuroplex_opentargets")
    print(f"  ✅ Appended {n_ot} GBA1 rows → neuroplex_opentargets")

print("\nDone. Verify counts:")
for tbl in ["neuroplex_gnomad", "neuroplex_uniprot", "neuroplex_kegg", "neuroplex_opentargets"]:
    cnt = spark.sql(f"SELECT COUNT(*) AS n FROM {CAT}.{SCH}.{tbl} WHERE gene_symbol = 'GBA1'").collect()[0]["n"]
    print(f"  {tbl:<30s}  GBA1 rows = {cnt}")

# COMMAND ----------

# DBTITLE 1,Deploy neuro-crispr-mcp
from databricks.sdk import WorkspaceClient
import time

w = WorkspaceClient()
app_name = "neuro-crispr-mcp"
source_path = "/Workspace/Users/andrew_forman@eisai.com/neuro-crispr-mcp"

# Step 1: check status before deploying
app = w.apps.get(app_name)
print(f"status        : {app.status}")
print(f"compute state : {app.compute_status.state if app.compute_status else 'N/A'}")
print(f"pending deploy: {app.pending_deployment}")

# Step 2: only deploy if RUNNING and no pending deployment
if str(app.status) == "AppStatus.RUNNING" and app.pending_deployment is None:
    print(f"\nDeploying from {source_path} ...")
    deployment = w.apps.deploy(
        app_name=app_name,
        source_code_path=source_path,
    )
    print(f"Deployment ID : {deployment.deployment_id}")
    print(f"Status        : {deployment.status}")
    print("Deployment started. Check app URL in ~2 min:")
    print(f"  https://neuro-crispr-mcp-3242147789106191.aws.databricksapps.com")
elif str(app.status) == "AppStatus.STOPPED":
    print("App is STOPPED — starting first ...")
    w.apps.start(app_name)
    print("Started. Wait ~2 min then re-run this cell to deploy.")
else:
    print(f"App state: {app.status} — check pending_deployment or wait for RUNNING.")

# COMMAND ----------

# DBTITLE 1,Deploy neuro-crispr-mcp
from databricks.sdk import WorkspaceClient

w        = WorkspaceClient()
APP_NAME = "neuro-crispr-mcp"
SRC_PATH = "/Workspace/Users/andrew_forman@eisai.com/neuro-crispr-mcp"

# Step 1: mandatory state check before any deploy/start
app = w.apps.get(APP_NAME)
status = str(app.status)
compute = str(app.compute_status.state) if app.compute_status else "N/A"
print(f"status        : {status}")
print(f"compute state : {compute}")
print(f"pending deploy: {app.pending_deployment}")

# Step 2: branch on state
if app.pending_deployment is not None:
    print("\u26a0\ufe0f  Deployment already in progress — poll apps.get() until clear, then re-run.")
elif "STOPPED" in status:
    print("App is STOPPED — starting first ...")
    w.apps.start(APP_NAME)
    print("  start() called. Wait ~2 min then re-run this cell.")
elif "RUNNING" in status:
    print(f"\nDeploying from {SRC_PATH} ...")
    d = w.apps.deploy(app_name=APP_NAME, source_code_path=SRC_PATH)
    print(f"  deployment_id : {d.deployment_id}")
    print(f"  deploy status : {d.status}")
    print(f"  App URL       : https://neuro-crispr-mcp-3242147789106191.aws.databricksapps.com")
else:
    print(f"Unexpected status {status!r} — inspect with w.apps.get() before retrying.")
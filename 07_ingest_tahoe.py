# Databricks notebook source
# /// script
# [tool.databricks.environment]
# environment_version = "5"
# ///
# DBTITLE 1,Install — huggingface_hub + datasets
# MAGIC %pip install -q huggingface_hub datasets pyarrow --upgrade

# COMMAND ----------

# DBTITLE 1,Config — environment, table targets, gene panel
import os, sys, json, time, logging
import pyarrow.parquet as pq
import pandas as pd
from pyspark.sql import functions as F
from pyspark.sql.types import (
    StructType, StructField, StringType, IntegerType, FloatType, DoubleType
)

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("tahoe_etl")

os.environ["NEUROPLEX_ENV"] = "val"
sys.path.insert(0, "/Workspace/Users/andrew_forman@eisai.com/neuro-crispr-mcp")
for mod in [k for k in sys.modules if k.startswith("config")]:
    del sys.modules[mod]
from config.neuroplex_config import load_config

CFG = load_config()
CAT, SCH = CFG.catalog, CFG.schema

TABLE_RAW       = f"{CAT}.{SCH}.tahoe_100m"
TABLE_CLUSTERED = f"{CAT}.{SCH}.tahoe_100m_clustered"
TABLE_GENE_VOCAB= f"{CAT}.{SCH}.tahoe_100m_gene_vocab"

# HuggingFace source
HF_REPO_ID  = "tahoebio/Tahoe-100M"
HF_REPO_TYPE = "dataset"

# NeuroPlex 31-gene panel (for gene-level filtering)
GENE_PANEL = [
    "PSEN1","PSEN2","APP","APOE","TREM2","BIN1",
    "CLU","CR1","PICALM","ABCA7","SORL1",
    "SNCA","LRRK2","PINK1","PRKN","PARK7","GBA1","VPS35",
    "SOD1","TARDBP","FUS","C9orf72","TBK1","GRN",
    "MAPT","HTT",
    "HCRT","HCRTR1","HCRTR2",
    "MECP2","FMR1",
]

print(f"✅ Environment : {CFG.environment}")
print(f"   Catalog     : {CAT}.{SCH}")
print(f"   Raw table   : {TABLE_RAW}")
print(f"   Clustered   : {TABLE_CLUSTERED}")
print(f"   Gene vocab  : {TABLE_GENE_VOCAB}")
print(f"   HF repo     : {HF_REPO_ID}")

# COMMAND ----------

# DBTITLE 1,Probe — list HuggingFace repo files and sizes
# Probe the tahoebio/Tahoe-100M repo to understand file layout and sizes before downloading.
# If the dataset is gated, set HUGGING_FACE_HUB_TOKEN in the environment first.

from huggingface_hub import list_repo_files, repo_info, hf_hub_url
import requests

# Check for HF token (needed if dataset is gated/private)
HF_TOKEN = os.environ.get("HUGGING_FACE_HUB_TOKEN", None)
if HF_TOKEN:
    print(f"  HF token   : ...{HF_TOKEN[-4:]} (authenticated)")
else:
    print("  HF token   : not set (attempting anonymous access)")

print(f"\nListing files in {HF_REPO_ID} ...\n")

try:
    files = list(
        list_repo_files(HF_REPO_ID, repo_type=HF_REPO_TYPE, token=HF_TOKEN)
    )
except Exception as exc:
    print(f"\u274c Failed to list files: {exc}")
    print("   If dataset is gated, set HUGGING_FACE_HUB_TOKEN env var.")
    files = []

# Group files by extension / prefix for readability
from collections import defaultdict
groups = defaultdict(list)
for f in files:
    ext = f.rsplit(".", 1)[-1] if "." in f else "other"
    groups[ext].append(f)

for ext, flist in sorted(groups.items()):
    print(f"  [{ext}] {len(flist)} files")
    for fn in flist[:8]:          # show first 8 per type
        print(f"      {fn}")
    if len(flist) > 8:
        print(f"      ... and {len(flist)-8} more")

print(f"\nTotal files : {len(files)}")

# COMMAND ----------

# DBTITLE 1,Schema probe — stream 1,000 rows to inspect columns and data types
# Read schema by downloading the FIRST parquet shard directly via hf_hub_download.
# Avoids the datasets streaming API (incompatible with upgraded huggingface_hub).
# Also fetches metadata/gene_vocabulary.json for the gene vocab table.

from huggingface_hub import hf_hub_download
import pyarrow.parquet as pq
import json as _json

# ---- 1. Identify the first parquet shard from the file listing ----
parquet_files = sorted([f for f in files if f.startswith("data/") and f.endswith(".parquet")])
print(f"Total parquet shards : {len(parquet_files)}")
print(f"First shard          : {parquet_files[0]}")
print(f"Last shard           : {parquet_files[-1]}")

# ---- 2. Download first shard and inspect schema + sample rows ----
print(f"\nDownloading {parquet_files[0]} ...")
first_shard_path = hf_hub_download(
    repo_id=HF_REPO_ID, filename=parquet_files[0],
    repo_type=HF_REPO_TYPE, token=HF_TOKEN,
)
pf = pq.ParquetFile(first_shard_path)
print(f"  Rows in shard       : {pf.metadata.num_rows:,}")
print(f"  Row groups          : {pf.metadata.num_row_groups}")
print(f"  Schema ({len(pf.schema_arrow)} fields):")
for field in pf.schema_arrow:
    print(f"    {field.name:<40s}  {str(field.type)}")

# Read 1,000 rows as pandas for value inspection
sample_df = pf.read_row_group(0).to_pandas().head(1000)
print(f"\nSample values (first 1,000 rows):")
for col in sample_df.columns:
    n_null = sample_df[col].isna().sum()
    sample_val = sample_df[col].dropna().iloc[0] if not sample_df[col].dropna().empty else None
    print(f"  {col:<40s}  dtype={str(sample_df[col].dtype):<14s}  "
          f"sample={repr(str(sample_val)[:60])}  nulls={n_null}")

# ---- 3. Fetch gene vocabulary ----
print("\n--- Gene Vocabulary ---")
gene_vocab_local = hf_hub_download(
    repo_id=HF_REPO_ID, filename="metadata/gene_vocabulary.json",
    repo_type=HF_REPO_TYPE, token=HF_TOKEN,
)
with open(gene_vocab_local) as fh:
    gene_vocab_raw = _json.load(fh)

# gene_vocabulary.json can be a list of names or a dict {gene: id}
if isinstance(gene_vocab_raw, list):
    gene_names = gene_vocab_raw
elif isinstance(gene_vocab_raw, dict):
    gene_names = list(gene_vocab_raw.keys())
else:
    gene_names = []
print(f"  Gene vocabulary size : {len(gene_names):,}")
print(f"  First 10 genes       : {gene_names[:10]}")
panel_in_vocab = [g for g in GENE_PANEL if g in set(gene_names)]
print(f"  NeuroPlex panel genes in vocab : {len(panel_in_vocab)}/31")
if len(panel_in_vocab) < 31:
    missing = [g for g in GENE_PANEL if g not in set(gene_names)]
    print(f"  Missing from vocab : {missing}")

# COMMAND ----------

# DBTITLE 1,Column mapping — normalise column names to server.py schema
# server.py expects: drug, moa, cell_line_id, cluster
# Tahoe-100M may use different column names — define the mapping here.
# UPDATE these after running the schema probe above.

# ---- ADJUST THESE BASED ON SCHEMA PROBE OUTPUT ----
COL_DRUG      = "drug"           # drug/compound/perturbation name
COL_MOA       = "moa-fine"        # Tahoe-100M column name (hyphenated — renamed to 'moa' on write)
COL_CELL_LINE = "cell_line_id"   # cell line identifier
COL_CLUSTER   = "plate"          # Tahoe-100M cluster proxy (no leiden column — plate used as cluster)
COL_DOSE      = "dose_um"        # dose in uM (optional — used in tahoe_100m raw table)
# ---- END ADJUSTABLE SECTION ----

# Auto-detect if sample_df is available
if 'sample_df' in dir():
    actual_cols = set(sample_df.columns)
    for attr_name, default_col in [
        ("COL_DRUG",      ["drug", "perturbation", "pert_iname", "compound"]),
        ("COL_MOA",       ["moa", "moa_id", "mechanism_of_action", "target"]),
        ("COL_CELL_LINE", ["cell_line_id", "cell_line", "cell_id", "line"]),
        ("COL_CLUSTER",   ["leiden", "cluster", "seurat_clusters", "cluster_id"]),
        ("COL_DOSE",      ["dose_um", "dose", "concentration", "pert_dose"]),
    ]:
        for candidate in default_col:
            if candidate in actual_cols:
                exec(f"{attr_name} = '{candidate}'")
                break

    print("Column mapping (auto-detected + overridable):")
    for attr in ["COL_DRUG", "COL_MOA", "COL_CELL_LINE", "COL_CLUSTER", "COL_DOSE"]:
        val = eval(attr)
        found = val in actual_cols
        print(f"  {attr:<16s} = {val!r:<30s}  {'\u2705' if found else '\u274c NOT FOUND in sample'}")
else:
    print("Run the schema probe cell first, then re-run this cell.")

# COMMAND ----------

# DBTITLE 1,Fetch drug annotations — MOA, target, approval status
# If MOA is NOT in the cell-level data, we need a separate drug annotation file.
# Common sources:
#  - A companion file in the HF repo (e.g. drug_metadata.csv / drug_annotations.parquet)
#  - PRISM/Broad drug repurposing hub annotations
#  - Or MOA may already be in the dataset (COL_MOA column present)

from huggingface_hub import hf_hub_download
import io

DRUG_ANNOT: dict[str, str] = {}   # drug_name → moa (populated below)

# ---- Step 1: Look for annotation file in the HF repo ----
ANNOT_CANDIDATES = [
    "drug_annotations.csv", "drug_annotations.parquet",
    "drug_metadata.csv",    "drug_metadata.parquet",
    "compounds.csv",        "compound_annotations.csv",
    "metadata/drugs.csv",   "data/drug_annotations.csv",
]
annot_file = None
if 'files' in dir():
    for candidate in ANNOT_CANDIDATES:
        if any(candidate in f for f in files):
            annot_file = next(f for f in files if candidate in f)
            print(f"Found annotation file: {annot_file}")
            break

if annot_file:
    try:
        local_path = hf_hub_download(
            repo_id=HF_REPO_ID, filename=annot_file,
            repo_type=HF_REPO_TYPE, token=HF_TOKEN
        )
        if annot_file.endswith(".parquet"):
            df_annot = pd.read_parquet(local_path)
        else:
            df_annot = pd.read_csv(local_path)
        print(f"  Loaded {len(df_annot):,} drug annotations")
        print(f"  Columns: {list(df_annot.columns)}")
        # Build drug → MOA mapping (try multiple column name conventions)
        drug_col = next((c for c in df_annot.columns if "drug" in c.lower() or "compound" in c.lower() or "name" in c.lower()), None)
        moa_col  = next((c for c in df_annot.columns if "moa" in c.lower() or "mechanism" in c.lower() or "target" in c.lower()), None)
        if drug_col and moa_col:
            DRUG_ANNOT = dict(zip(df_annot[drug_col].str.lower(), df_annot[moa_col].fillna("")))
            print(f"  MOA mapping: {len(DRUG_ANNOT):,} entries  (drug_col={drug_col!r}, moa_col={moa_col!r})")
            for k, v in list(DRUG_ANNOT.items())[:5]:
                print(f"    {k!r:30s} → {v!r}")
    except Exception as exc:
        print(f"  ⚠️  Could not load annotation file: {exc}")
else:
    print("⚠️  No drug annotation file found in repo.")
    print("   If moa is already a column in the dataset (check schema probe), DRUG_ANNOT stays empty.")
    print("   If needed, MOA can be sourced from PRISM repurposing hub:")
    print("     https://www.broadinstitute.org/drug-repurposing-hub")

print(f"\nDRUG_ANNOT entries: {len(DRUG_ANNOT):,}")

# COMMAND ----------

# DBTITLE 1,ETL — download all shards via hf_hub_download + pyarrow, write tahoe_100m
# Tahoe-100M incremental ETL.
# Downloads each parquet shard via hf_hub_download, reads only 8 metadata columns
# (skipping the large sparse genes/expressions lists), then appends to Delta.
# Batches 50 shards per Delta write. Resume-safe: estimates progress from row count.
#
# Runtime: ~30-120 min depending on shard file sizes and network speed.

import os, sys, time, warnings
import pyarrow.parquet as pq
import pandas as pd
from huggingface_hub import hf_hub_download
from pyspark.sql.types import StructType, StructField, StringType

warnings.filterwarnings("ignore")

for mod in [k for k in sys.modules if k.startswith("config")]:
    del sys.modules[mod]
from config.neuroplex_config import load_config
_cfg = load_config()

TABLE_META  = f"{_cfg.catalog}.{_cfg.schema}.tahoe_100m"
BATCH_SIZE  = 50     # shards per Delta write
MAX_SHARDS  = None   # None = all 3388; set an integer to cap for a pilot run

# ---- Shard list from cell 3 probe ----
if "files" not in dir():
    raise RuntimeError("Run cell 3 (Probe) first to populate `files`.")

parquet_shards = sorted([
    f for f in files
    if f.endswith(".parquet")
    and not any(x in f.lower() for x in
                ["var", "obs", "vocab", "gene", "drug_annot", "annot", "metadata"])
])
if not parquet_shards:              # fallback: any parquet
    parquet_shards = sorted([f for f in files if f.endswith(".parquet")])

print(f"Shard files found : {len(parquet_shards):,}")

if MAX_SHARDS:
    parquet_shards = parquet_shards[:MAX_SHARDS]
    print(f"Capped at         : {MAX_SHARDS} shards")

# ---- Resume: skip already-ingested batches ----
try:
    existing_rows = spark.sql(f"SELECT COUNT(*) AS n FROM {TABLE_META}").collect()[0]["n"]
    est_shards    = int(existing_rows / 28225)           # ~28K rows per shard
    skip_shards   = (est_shards // BATCH_SIZE) * BATCH_SIZE
    if skip_shards > 0:
        print(f"Existing rows     : {existing_rows:,}  (~{est_shards} shards done)")
        print(f"Resuming from shard {skip_shards + 1}")
        parquet_shards = parquet_shards[skip_shards:]
    else:
        # Placeholder table exists with wrong schema — drop it so first write creates fresh
        print("Dropping placeholder table (schema mismatch) and starting fresh ETL ...")
        spark.sql(f"DROP TABLE IF EXISTS {TABLE_META}")
        existing_rows = 0
except Exception:
    existing_rows = 0
    skip_shards   = 0
    print("Starting fresh ETL")

# ---- Column spec (exclude genes / expressions list columns) ----
READ_COLS = [
    COL_DRUG, COL_MOA, COL_CELL_LINE, COL_CLUSTER,
    "canonical_smiles", "pubchem_cid", "sample", "BARCODE_SUB_LIB_ID",
]
RENAME_MAP = {
    COL_MOA:               "moa",
    COL_CLUSTER:           "cluster",
    COL_DRUG:              "drug",
    COL_CELL_LINE:         "cell_line_id",
    "canonical_smiles":    "canonical_smiles",
    "pubchem_cid":         "pubchem_cid",
    "sample":              "sample",
    "BARCODE_SUB_LIB_ID": "barcode",
}
SCHEMA_COLS = ["drug", "moa", "cell_line_id", "cluster",
               "canonical_smiles", "pubchem_cid", "sample", "barcode"]

SPARK_SCHEMA = StructType([
    StructField(c, StringType(), True) for c in SCHEMA_COLS
])


def _read_shard(fname: str) -> pd.DataFrame | None:
    """Download shard to /tmp, read metadata cols only, delete file."""
    tmp = None
    try:
        tmp = hf_hub_download(
            repo_id=HF_REPO_ID, filename=fname,
            repo_type=HF_REPO_TYPE, token=HF_TOKEN,
        )
        pf        = pq.ParquetFile(tmp)
        available = [c for c in READ_COLS if c in pf.schema_arrow.names]
        df        = pf.read(columns=available).to_pandas()
        df        = df.rename(columns={k: v for k, v in RENAME_MAP.items() if k in df.columns})
        for col in SCHEMA_COLS:
            if col not in df.columns:
                df[col] = None
        return df[SCHEMA_COLS]
    except Exception as exc:
        print(f"    ⚠️  {fname}: {exc}")
        return None
    finally:
        if tmp and os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass


# ---- Main ingestion loop ----
t0         = time.time()
total_rows = existing_rows
n_shards   = len(parquet_shards)
n_batches  = (n_shards + BATCH_SIZE - 1) // BATCH_SIZE
base_batch = skip_shards // BATCH_SIZE

print(f"\nIngesting {n_shards:,} shards in {n_batches} batches of {BATCH_SIZE} ...")

for b in range(n_batches):
    batch_files = parquet_shards[b * BATCH_SIZE : (b + 1) * BATCH_SIZE]
    b_t0 = time.time()

    dfs = [_read_shard(f) for f in batch_files]
    dfs = [d for d in dfs if d is not None and not d.empty]
    if not dfs:
        continue

    batch_pd = pd.concat(dfs, ignore_index=True).astype(str)
    batch_pd = batch_pd.where(batch_pd != "None", other=None)
    batch_pd = batch_pd.where(batch_pd != "nan",  other=None)

    write_mode = "overwrite" if (b == 0 and skip_shards == 0) else "append"
    w = spark.createDataFrame(batch_pd, schema=SPARK_SCHEMA).write.format("delta")
    if write_mode == "overwrite":
        w = w.option("overwriteSchema", "true")
    w.mode(write_mode).saveAsTable(TABLE_META)

    total_rows  += len(batch_pd)
    b_elapsed    = time.time() - b_t0
    remaining    = n_batches - b - 1
    eta_min      = remaining * b_elapsed / 60

    print(
        f"  Batch {base_batch + b + 1:3d}/{base_batch + n_batches}: "
        f"+{len(batch_pd):>6,} | total {total_rows:>10,} | "
        f"{b_elapsed:.1f}s | ETA ~{eta_min:.0f} min"
    )

print(f"\n✅ tahoe_100m: {total_rows:,} rows in {(time.time()-t0)/60:.1f} min")

# COMMAND ----------

# DBTITLE 1,Build tahoe_100m_clustered via CTAS + GRANT SELECT to app SP
# Build the table all three server.py Tahoe tools query.
# CTAS is fast (server-side operation) and avoids re-downloading any data.

import os, sys
os.environ["NEUROPLEX_ENV"] = "val"
for mod in [k for k in sys.modules if k.startswith("config")]:
    del sys.modules[mod]
from config.neuroplex_config import load_config
_cfg = load_config()

TABLE_META  = f"{_cfg.catalog}.{_cfg.schema}.tahoe_100m"
TABLE_CLU   = f"{_cfg.catalog}.{_cfg.schema}.tahoe_100m_clustered"
SP_CLIENT_ID = "79e31bca-0a8a-4fa8-b858-976114fb1f9b"

print("Building tahoe_100m_clustered via CTAS ...")
spark.sql(f"""
    CREATE OR REPLACE TABLE {TABLE_CLU}
    AS
    SELECT drug, moa, cell_line_id, cluster
    FROM   {TABLE_META}
    WHERE  drug IS NOT NULL
""")
n = spark.sql(f"SELECT COUNT(*) AS n FROM {TABLE_CLU}").collect()[0]["n"]
print(f"  {n:,} rows → tahoe_100m_clustered")

# Optimize for GROUP BY queries (drug, moa, cell_line_id)
print("Optimising + Z-ordering ...")
spark.sql(f"OPTIMIZE {TABLE_CLU} ZORDER BY (drug, moa, cell_line_id)")
print("  Done")

# Grant SELECT to app SP
try:
    for tbl in [TABLE_META, TABLE_CLU]:
        spark.sql(f"GRANT SELECT ON TABLE {tbl} TO `{SP_CLIENT_ID}`")
        print(f"  GRANT SELECT → {tbl.split('.')[-1]}")
except Exception as _grant_err:
    print(f"  ⚠️  GRANT skipped (permission error): {_grant_err}")

print("✅ tahoe_100m_clustered ready")

# COMMAND ----------

# DBTITLE 1,Fetch gene vocabulary from HuggingFace var file
import pandas as pd

# Gene vocabulary: the list of genes measured in Tahoe-100M.
# Typically stored in a var.parquet / gene_names.txt / features.tsv.gz in the HF repo.

GENE_VOCAB_CANDIDATES = [
    "var.parquet", "var.csv",
    "genes.parquet", "genes.csv",
    "gene_names.txt", "gene_names.csv",
    "features.tsv.gz", "features.tsv",
    "data/var.parquet", "metadata/genes.csv",
]

vocab_file = None
if 'files' in dir():
    for candidate in GENE_VOCAB_CANDIDATES:
        if any(candidate in f for f in files):
            vocab_file = next(f for f in files if candidate in f)
            print(f"Found gene vocab file: {vocab_file}")
            break

gene_vocab_df = None
if vocab_file:
    try:
        local_path = hf_hub_download(
            repo_id=HF_REPO_ID, filename=vocab_file,
            repo_type=HF_REPO_TYPE, token=HF_TOKEN
        )
        if vocab_file.endswith(".parquet"):
            gene_vocab_df = pd.read_parquet(local_path)
        elif vocab_file.endswith((".tsv.gz", ".tsv")):
            gene_vocab_df = pd.read_csv(local_path, sep="\t", header=None, names=["gene_id", "gene_name", "feature_type"])
        else:
            gene_vocab_df = pd.read_csv(local_path)
        print(f"  {len(gene_vocab_df):,} genes  |  columns: {list(gene_vocab_df.columns)}")
        print(gene_vocab_df.head(5).to_string())
    except Exception as exc:
        print(f"  ⚠️  Could not load vocab file: {exc}")
else:
    print("⚠️  No gene vocab file found in HF repo.")
    print("   Will infer gene list from dataset column names (if expression columns exist).")
    # Fallback: use gene panel + mark as subset
    gene_vocab_df = pd.DataFrame({"gene_name": GENE_PANEL, "source": "neuroplex_panel"})

# COMMAND ----------

# DBTITLE 1,Write tahoe_100m_gene_vocab to Delta
# Normalise and write tahoe_100m_gene_vocab.
# Required columns: gene_name (STRING).  Optional: gene_id, ensembl_id, is_panel_gene.

if gene_vocab_df is None:
    print("  No gene vocabulary available — skipping.")
else:
    vocab = gene_vocab_df.copy()

    # Normalise column names
    rename_map = {}
    for col in vocab.columns:
        lc = col.lower()
        if any(k in lc for k in ["gene_name", "gene_symbol", "symbol", "name"]):
            rename_map[col] = "gene_name"
        elif any(k in lc for k in ["gene_id", "ensembl", "ensg"]):
            rename_map[col] = "ensembl_id"
    vocab = vocab.rename(columns=rename_map)

    if "gene_name" not in vocab.columns:
        # Use first column as gene_name
        vocab = vocab.rename(columns={vocab.columns[0]: "gene_name"})

    # Add flag: is this gene in the NeuroPlex panel?
    vocab["is_panel_gene"] = vocab["gene_name"].isin(GENE_PANEL)
    panel_count = vocab["is_panel_gene"].sum()

    print(f"  {len(vocab):,} genes  |  {panel_count} in NeuroPlex panel")
    print(f"  Columns: {list(vocab.columns)}")

    df_vocab = spark.createDataFrame(vocab.astype(str))
    df_vocab.write.format("delta").mode("overwrite").option("overwriteSchema", "true").saveAsTable(TABLE_GENE_VOCAB)
    print(f"  ✅ Written {df_vocab.count():,} rows → {TABLE_GENE_VOCAB}")

# COMMAND ----------

# DBTITLE 1,Validation — top drugs by cell count
# MAGIC %sql
# MAGIC -- Mirrors list_tahoe_drugs query in server.py
# MAGIC -- Expect: ~500-1000 drugs, largest MOA classes at top
# MAGIC SELECT
# MAGIC     drug,
# MAGIC     moa,
# MAGIC     COUNT(*)                    AS n_cells,
# MAGIC     COUNT(DISTINCT cell_line_id) AS n_cell_lines
# MAGIC FROM dhbl_discovery_us_val.genesis_schema.tahoe_100m_clustered
# MAGIC WHERE drug IS NOT NULL
# MAGIC GROUP BY drug, moa
# MAGIC ORDER BY n_cells DESC
# MAGIC LIMIT 20

# COMMAND ----------

# DBTITLE 1,Validation — MOA distribution
# MAGIC %sql
# MAGIC -- Mirrors query_tahoe_moa_distribution in server.py
# MAGIC SELECT
# MAGIC     moa,
# MAGIC     COUNT(DISTINCT drug)         AS n_drugs,
# MAGIC     COUNT(*)                     AS n_cells,
# MAGIC     COUNT(DISTINCT cell_line_id) AS n_cell_lines
# MAGIC FROM dhbl_discovery_us_val.genesis_schema.tahoe_100m_clustered
# MAGIC WHERE moa IS NOT NULL AND moa != 'unknown'
# MAGIC GROUP BY moa
# MAGIC ORDER BY n_drugs DESC
# MAGIC LIMIT 20

# COMMAND ----------

# DBTITLE 1,Validation — row counts across all three Tahoe tables
# MAGIC %sql
# MAGIC SELECT 'tahoe_100m'           AS tbl, COUNT(*) AS rows FROM dhbl_discovery_us_val.genesis_schema.tahoe_100m
# MAGIC UNION ALL
# MAGIC SELECT 'tahoe_100m_clustered' AS tbl, COUNT(*) AS rows FROM dhbl_discovery_us_val.genesis_schema.tahoe_100m_clustered
# MAGIC UNION ALL
# MAGIC SELECT 'tahoe_100m_gene_vocab'AS tbl, COUNT(*) AS rows FROM dhbl_discovery_us_val.genesis_schema.tahoe_100m_gene_vocab
# MAGIC ORDER BY tbl

# COMMAND ----------

# DBTITLE 1,Probe — perturbai/wholebrain_crispr_atlas schemas
# Memory-safe schema probe for perturbai/wholebrain_crispr_atlas.
# Reads only schema + first row group (NOT full files) to avoid OOM.

import os, pyarrow.parquet as pq
import pandas as pd
from huggingface_hub import list_repo_files, hf_hub_download
from collections import defaultdict

CRISPR_REPO = "perturbai/wholebrain_crispr_atlas"
HF_TOKEN    = os.environ.get("HUGGING_FACE_HUB_TOKEN", None)

files_crispr = list(list_repo_files(CRISPR_REPO, repo_type="dataset", token=HF_TOKEN))
data_shards  = sorted([f for f in files_crispr if f.startswith("data/") and f.endswith(".parquet")])
de_file      = next(f for f in files_crispr if "wilcoxon_de_results" in f)
ndd_csvs     = [f for f in files_crispr if "ndd_perturbation" in f and f.split(".")[-1] in ("csv","tsv")]
print(f"Data shards: {len(data_shards)}   DE file: {de_file.split('/')[-1]}")
print(f"NDD files  : {[f.split('/')[-1] for f in ndd_csvs]}")

# ---- Cell data schema (first shard, first 5 rows) ----
print("\n=== Cell data shard ===")
p = hf_hub_download(CRISPR_REPO, data_shards[0], repo_type="dataset", token=HF_TOKEN)
pf = pq.ParquetFile(p)
print(f"  Rows/shard: {pf.metadata.num_rows:,}  |  Fields: {len(pf.schema_arrow)}")
for field in pf.schema_arrow:
    print(f"    {field.name:<40s}  {str(field.type)}")
small = pf.read_row_group(0, columns=list(pf.schema_arrow.names)[:10]).to_pandas().head(5)
print(small.to_string())
del pf, small

# ---- DE results: schema + 5 rows ONLY, no full load ----
print("\n=== DE results ===")
p2 = hf_hub_download(CRISPR_REPO, de_file, repo_type="dataset", token=HF_TOKEN)
pf2 = pq.ParquetFile(p2)
print(f"  Row groups: {pf2.metadata.num_row_groups}  |  Total rows: {pf2.metadata.num_rows:,}")
print(f"  Fields    : {list(pf2.schema_arrow.names)}")
de_sample = pf2.read_row_group(0).to_pandas().head(5)
print(de_sample.to_string())
del pf2, de_sample

# ---- NDD CSV schemas (small files) ----
print("\n=== NDD gene files ===")
for csv_f in ndd_csvs[:4]:
    try:
        lp = hf_hub_download(CRISPR_REPO, csv_f, repo_type="dataset", token=HF_TOKEN)
        sep = "\t" if csv_f.endswith(".tsv") else ","
        df_n = pd.read_csv(lp, sep=sep, nrows=200)
        print(f"\n  {csv_f.split('/')[-1]:55s}  ~{len(df_n)} rows  cols={list(df_n.columns)}")
        print(df_n.head(3).to_string())
    except Exception as e:
        print(f"  ERROR {csv_f}: {e}")
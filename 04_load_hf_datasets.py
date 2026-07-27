# Databricks notebook source
# /// script
# [tool.databricks.environment]
# environment_version = "5"
# ///
# DBTITLE 1,HuggingFace Dataset Loader — Tahoe-100M + CRISPR Atlas
# MAGIC %md
# MAGIC # HuggingFace → Delta ETL: Tahoe-100M + CRISPR Atlas
# MAGIC
# MAGIC | Dataset | HuggingFace repo | Tables written |
# MAGIC |---|---|---|
# MAGIC | CRISPR Atlas | `perturbai/wholebrain_crispr_atlas` | `crispr_atlas_gene_metadata`, `crispr_atlas_cell_metadata`, `crispr_atlas_ndd_genes`, `crispr_atlas_diff_expr`, `wholebrain_crispr_atlas` |
# MAGIC | Tahoe drug DE | `tahoebio/tahoe-de-rhaister` | `tahoe_100m_clustered` |
# MAGIC | Tahoe metadata | `tahoebio/Tahoe-100M` | `tahoe_100m_gene_vocab`, `tahoe_100m` (drug catalog) |
# MAGIC
# MAGIC Run cells top-to-bottom. Large file cells (Atlas DE, Tahoe streaming) are safe to re-run — they append/merge.

# COMMAND ----------

# DBTITLE 1,Setup — install deps, configure env
# MAGIC %pip install -q huggingface_hub pyarrow
# MAGIC
# MAGIC import os, sys, time, json, logging
# MAGIC os.environ.setdefault("NEUROPLEX_ENV", "prod")
# MAGIC
# MAGIC for mod in list(sys.modules):
# MAGIC     if "ingestion" in mod or "neuroplex_config" in mod:
# MAGIC         del sys.modules[mod]
# MAGIC
# MAGIC sys.path.insert(0, "/Workspace/Users/andrew_forman@eisai.com/neuro-crispr-mcp")
# MAGIC
# MAGIC from ingestion.source_registry import TARGET_CATALOG, TARGET_SCHEMA
# MAGIC
# MAGIC CAT    = TARGET_CATALOG
# MAGIC SCH    = TARGET_SCHEMA
# MAGIC PREFIX = f"{CAT}.{TARGET_SCHEMA}"
# MAGIC
# MAGIC ATLAS_REPO  = "perturbai/wholebrain_crispr_atlas"
# MAGIC TAHOE_REPO  = "tahoebio/Tahoe-100M"
# MAGIC TAHOE_DE_REPO = "tahoebio/tahoe-de-rhaister"
# MAGIC
# MAGIC print(f"Target catalog : {CAT}")
# MAGIC print(f"Target schema  : {SCH}")
# MAGIC print(f"Full prefix    : {PREFIX}")

# COMMAND ----------

# DBTITLE 1,Atlas — crispr_atlas_gene_metadata (with fallback to DE extraction)
from huggingface_hub import hf_hub_download
import pandas as pd, pyarrow.parquet as pq

TABLE = f"{PREFIX}.crispr_atlas_gene_metadata"
gene_df = None

# Try known metadata paths first; the repo may or may not have a metadata/ dir
for candidate in [
    "metadata/gene_metadata.parquet",
    "metadata/var.parquet",
    "analysis/genes_of_interest_blomen.csv",
]:
    try:
        p = hf_hub_download(ATLAS_REPO, candidate, repo_type="dataset")
        gene_df = pd.read_parquet(p) if candidate.endswith(".parquet") else pd.read_csv(p)
        print(f"✅ Loaded from {candidate}: {gene_df.shape}")
        print(f"   Columns: {list(gene_df.columns)}")
        print(gene_df.head(3).to_string())
        break
    except Exception as e:
        print(f"   {candidate}: {e}")

if gene_df is None:
    # metadata/ does not exist — derive gene list from the Wilcoxon DE results
    # (always present at the confirmed path)
    print("\nFalling back: extracting gene metadata from DE results…")
    DE_FILE = (
        "analysis/2603_shi_manuscript/diff_expr_genes/"
        "wilcoxon_de_results_group_name_final_data_no_multi_guide_cells.parquet"
    )
    p_de = hf_hub_download(ATLAS_REPO, DE_FILE, repo_type="dataset")

    # Peek at schema before loading the full file
    schema_de = pq.read_schema(p_de)
    all_cols = [f.name for f in schema_de]
    print(f"   DE columns ({len(all_cols)}): {all_cols}")

    df_de = pd.read_parquet(p_de)
    print(f"   DE shape: {df_de.shape}")
    print(f"   Sample:\n{df_de.head(3).to_string()}")

    # Scanpy convention: 'names' = gene symbol.  Fall back to first column.
    gene_col = next(
        (c for c in df_de.columns
         if c in ("names", "gene", "gene_symbol", "perturbation", "gene_target")),
        df_de.columns[0],
    )
    print(f"\n   Using column '{gene_col}' as gene identifier")

    # Unique gene targets → gene_metadata table
    unique_genes = df_de[gene_col].dropna().unique()
    gene_df = pd.DataFrame({
        "gene_symbol": unique_genes,
        "source":      "derived_from_wilcoxon_de",
    })
    print(f"   Extracted {len(gene_df):,} unique gene targets")

sdf = spark.createDataFrame(gene_df.astype(str))
sdf.write.mode("overwrite").option("overwriteSchema", "true").saveAsTable(TABLE)
cnt = spark.sql(f"SELECT COUNT(*) AS n FROM {TABLE}").collect()[0]["n"]
print(f"\n✅ {TABLE}: {cnt:,} rows")

# COMMAND ----------

# DBTITLE 1,Atlas — crispr_atlas_ndd_genes (NDD gene lists)
import io, requests

TABLE = f"{PREFIX}.crispr_atlas_ndd_genes"

ndd_files = [
    "analysis/2603_shi_manuscript/ndd_perturbation/unique_genes_all.csv",
    "analysis/2603_shi_manuscript/ndd_perturbation/ASD_rare.csv",
    "analysis/2603_shi_manuscript/ndd_perturbation/GWS_results.csv",
]

dfs = []
for fname in ndd_files:
    try:
        p = hf_hub_download(ATLAS_REPO, fname, repo_type="dataset")
        df = pd.read_csv(p)
        df["source_file"] = fname.split("/")[-1]
        dfs.append(df)
        print(f"  ✅ {fname.split('/')[-1]}: {df.shape}")
    except Exception as e:
        print(f"  ⚠️  {fname.split('/')[-1]}: {e}")

if dfs:
    combined = pd.concat(dfs, ignore_index=True)
    sdf = spark.createDataFrame(combined.astype(str))  # all strings — schemas differ per file
    sdf.write.mode("overwrite").option("overwriteSchema", "true").saveAsTable(TABLE)
    cnt = spark.sql(f"SELECT COUNT(*) AS n FROM {TABLE}").collect()[0]["n"]
    print(f"\n✅ {TABLE}: {cnt:,} rows")

# COMMAND ----------

# DBTITLE 1,Atlas — crispr_atlas_cell_metadata (fallback to WB8588 stream)
from huggingface_hub import hf_hub_download, list_repo_files
import pandas as pd, pyarrow.parquet as pq

TABLE = f"{PREFIX}.crispr_atlas_cell_metadata"
done = False

# Try known metadata paths first
for candidate in ["metadata/all_obs.parquet", "metadata/obs.parquet"]:
    try:
        print(f"Trying {candidate} (~470 MB expected)…")
        t0 = time.time()
        p = hf_hub_download(ATLAS_REPO, candidate, repo_type="dataset")
        print(f"  Downloaded in {time.time()-t0:.1f}s")
        # Serverless blocks file:// reads — use pyarrow batch iteration instead
        reader = pq.ParquetFile(p)
        n_rows = reader.metadata.num_rows
        print(f"  {n_rows:,} rows — reading in 250K-row batches…")
        for i, batch in enumerate(reader.iter_batches(batch_size=250_000)):
            df_b = batch.to_pandas()
            for c in df_b.select_dtypes(include="category").columns:
                df_b[c] = df_b[c].astype(str)
            spark.createDataFrame(df_b).write \
                .mode("overwrite" if i == 0 else "append") \
                .option("overwriteSchema", "true") \
                .saveAsTable(TABLE)
            print(f"  batch {i+1}: {len(df_b):,} rows  ({(i+1)*250_000/n_rows*100:.0f}%)")
        done = True
        break
    except Exception as e:
        print(f"  {candidate}: {e}")

if not done:
    # metadata/ absent — stream from WB8588_* parquet shards
    print("\nmetadata/ not found — streaming WB8588 cell data files…")
    wb_files = sorted([
        f for f in list_repo_files(ATLAS_REPO, repo_type="dataset")
        if f.startswith("data/WB8588_") and f.endswith(".parquet")
    ])
    print(f"Found {len(wb_files)} WB8588 files")

    # Probe first file to discover schema and metadata columns
    p0 = hf_hub_download(ATLAS_REPO, wb_files[0], repo_type="dataset")
    schema0 = pq.read_schema(p0)
    all_cols = [f.name for f in schema0]
    print(f"\nWB8588 schema ({len(all_cols)} cols): {all_cols}")
    sample = pd.read_parquet(p0)
    print(f"Sample row:\n{sample.iloc[0].to_dict()}")

    import numpy as np
    meta_cols = [
        c for c in all_cols
        if not pd.api.types.is_float_dtype(sample[c]) or sample[c].nunique() < 200
    ]
    meta_keywords = {"cell", "type", "cluster", "umap", "barcode", "perturb",
                     "batch", "sample", "gene_target", "group", "leiden", "louvain"}
    explicit_meta = [
        c for c in all_cols
        if any(k in c.lower() for k in meta_keywords)
    ]
    keep_cols = list(set(meta_cols) | set(explicit_meta))
    print(f"\nKeeping {len(keep_cols)} metadata columns: {keep_cols}")

    # Stream all WB8588 files using pandas (serverless blocks file:// reads)
    for i, fname in enumerate(wb_files):
        try:
            p = hf_hub_download(ATLAS_REPO, fname, repo_type="dataset")
            df_shard = pd.read_parquet(p)
            shard_keep = [c for c in keep_cols if c in df_shard.columns]
            df_meta = df_shard[shard_keep] if shard_keep else df_shard
            for c in df_meta.select_dtypes(include="category").columns:
                df_meta[c] = df_meta[c].astype(str)
            spark.createDataFrame(df_meta).write \
                .mode("overwrite" if i == 0 else "append") \
                .option("overwriteSchema", "true") \
                .option("mergeSchema",     "true") \
                .saveAsTable(TABLE)
            if i % 20 == 0 or i == len(wb_files) - 1:
                interim = spark.sql(f"SELECT COUNT(*) AS n FROM {TABLE}").collect()[0]["n"]
                print(f"  [{i+1}/{len(wb_files)}] {fname.split('/')[-1]} — table: {interim:,} rows")
        except Exception as e:
            print(f"  ⚠️  {fname}: {e}")

cnt = spark.sql(f"SELECT COUNT(*) AS n FROM {TABLE}").collect()[0]["n"]
print(f"\n✅ {TABLE}: {cnt:,} rows")

# COMMAND ----------

# DBTITLE 1,Atlas — crispr_atlas_diff_expr + wholebrain_crispr_atlas (Wilcoxon DE, Spark only)
# 745M rows / 10.8 GB — too large for pandas or spark.read.parquet(file://).
# Strategy: pyarrow batch reader, 5 row-groups at a time (~5M rows),
# filter to pvals_adj < 0.05, write each batch to Delta in append mode.
# Mirrors result to wholebrain_crispr_atlas via CTAS at the end.
import pyarrow as pa, pyarrow.parquet as pq, gc

DE_FILE = (
    "analysis/2603_shi_manuscript/diff_expr_genes/"
    "wilcoxon_de_results_group_name_final_data_no_multi_guide_cells.parquet"
)
TABLE_DE   = f"{PREFIX}.crispr_atlas_diff_expr"
TABLE_MAIN = f"{PREFIX}.wholebrain_crispr_atlas"

RENAME = {
    "names": "de_gene", "logfoldchanges": "log2fc",
    "pvals": "pvalue",  "pvals_adj": "pvalue_adj",
    "n_pert_matched": "n_cells_ko", "n_ctrl_matched": "n_cells_ctrl",
}
PVAL_CUT = 0.05
BATCH_G  = 5   # row-groups per batch → ~5.25 M raw rows, ~420 MB in-memory peak

print("Locating Wilcoxon DE parquet (cached after earlier download)…")
local_de  = hf_hub_download(ATLAS_REPO, DE_FILE, repo_type="dataset")
reader    = pq.ParquetFile(local_de)
n_groups  = reader.metadata.num_rows  # note: num_row_groups for iteration
n_groups  = reader.metadata.num_row_groups
print(f"File: {reader.metadata.num_rows:,} rows | {n_groups} row groups | pval cutoff={PVAL_CUT}")

total_written = 0
t0 = time.time()

for batch_start in range(0, n_groups, BATCH_G):
    batch_end = min(batch_start + BATCH_G, n_groups)
    tables = [reader.read_row_group(g) for g in range(batch_start, batch_end)]
    df = pa.concat_tables(tables).to_pandas()

    df_sig = df[df["pvals_adj"] < PVAL_CUT].rename(columns=RENAME).copy()
    del df; gc.collect()

    if len(df_sig) > 0:
        mode = "overwrite" if total_written == 0 else "append"
        spark.createDataFrame(df_sig).write \
            .mode(mode).option("overwriteSchema", "true").saveAsTable(TABLE_DE)
        total_written += len(df_sig)
    del df_sig; gc.collect()

    elapsed = time.time() - t0
    print(f"  [{batch_end:>3}/{n_groups}] {batch_end/n_groups*100:>5.1f}%  "
          f"written: {total_written:>10,}  {elapsed:.0f}s")

# Mirror to wholebrain_crispr_atlas via CTAS (single scan, no re-download)
print("\nMirroring to wholebrain_crispr_atlas via CTAS…")
spark.sql(f"CREATE OR REPLACE TABLE {TABLE_MAIN} AS SELECT * FROM {TABLE_DE}")

c1 = spark.sql(f"SELECT COUNT(*) AS n FROM {TABLE_DE}").collect()[0]["n"]
c2 = spark.sql(f"SELECT COUNT(*) AS n FROM {TABLE_MAIN}").collect()[0]["n"]
print(f"\n✅ crispr_atlas_diff_expr  : {c1:,} rows")
print(f"✅ wholebrain_crispr_atlas : {c2:,} rows")
print(f"Total time: {time.time()-t0:.0f}s")

# COMMAND ----------

# DBTITLE 1,Tahoe — tahoe_100m_gene_vocab (Entrez ID → gene symbol)
import numpy as np

TABLE = f"{PREFIX}.tahoe_100m_gene_vocab"

# The Tahoe-100M stores expression as parallel arrays: genes[i] = Entrez ID, expressions[i] = z-score.
# Extract the gene vocabulary from the first shard (all shards share the same gene set).
print("Downloading Tahoe shard 0 for gene vocab…")
path = hf_hub_download(TAHOE_REPO, "data/train-00000-of-03388.parquet", repo_type="dataset")
df   = pd.read_parquet(path, columns=["genes"])

# Each row's genes array = list/ndarray of Entrez IDs present in that cell (non-zero)
# Collect the union of all gene IDs across all cells in this shard
all_genes = set()
for arr in df["genes"]:
    all_genes.update(arr.tolist() if hasattr(arr, 'tolist') else list(arr))

print(f"Unique gene IDs in shard 0: {len(all_genes):,}")

# Create vocab table: gene_token_id (Entrez int) → will resolve to symbol via NCBI
vocab_df = pd.DataFrame({"gene_id": sorted(all_genes)}).astype({"gene_id": "int64"})
vocab_df["gene_symbol"] = None   # placeholder; enrich with NCBI lookup if needed
vocab_df["source"] = "tahoe_100m_shard0"

print(f"Vocab size: {len(vocab_df):,}")
print(vocab_df.head(5))

sdf = spark.createDataFrame(vocab_df)
sdf.write.mode("overwrite").option("overwriteSchema", "true").saveAsTable(TABLE)
cnt = spark.sql(f"SELECT COUNT(*) AS n FROM {TABLE}").collect()[0]["n"]
print(f"\n✅ {TABLE}: {cnt:,} rows")

# COMMAND ----------

# DBTITLE 1,Tahoe — tahoe_100m drug catalog (sampled metadata sweep, fast)
from huggingface_hub import list_repo_files, hf_hub_download
import pandas as pd

TABLE = f"{PREFIX}.tahoe_100m_clustered"

# tahoe-de-rhaister: cell_eval/plate_plate*.parquet
# Wide format: [cell_line, treatment, GENE1, GENE2, … ~2000 gene cols]
# Melt → long: cell_line, drug, gene_symbol, de_score
# Keep top 200 DE genes per treatment×cell_line to control table size.

plate_files = sorted([
    f for f in list_repo_files(TAHOE_DE_REPO, repo_type="dataset")
    if f.startswith("cell_eval/plate_plate") and f.endswith(".parquet")
])
print(f"Found {len(plate_files)} plate files in tahoe-de-rhaister")

total_rows = 0
for i, fname in enumerate(plate_files):
    print(f"  [{i+1}/{len(plate_files)}] {fname}…", end=" ", flush=True)
    t0 = time.time()
    try:
        p = hf_hub_download(TAHOE_DE_REPO, fname, repo_type="dataset")
        df_wide = pd.read_parquet(p)

        meta_cols = ["cell_line", "treatment"]
        gene_cols = [c for c in df_wide.columns if c not in meta_cols]

        df_long = df_wide.melt(
            id_vars=meta_cols,
            value_vars=gene_cols,
            var_name="gene_symbol",
            value_name="de_score",
        ).dropna(subset=["de_score"])
        df_long["plate"] = fname.split("plate_")[-1].replace(".parquet", "")
        df_long["drug"]  = df_long["treatment"]

        df_long["abs_de"] = df_long["de_score"].abs()
        df_top = (
            df_long
            .sort_values("abs_de", ascending=False)
            .groupby(["cell_line", "drug"], group_keys=False)
            .head(200)
            .drop(columns=["abs_de"])
        )

        sdf = spark.createDataFrame(df_top)
        mode = "overwrite" if i == 0 else "append"
        sdf.write.mode(mode).option("overwriteSchema", "true").saveAsTable(TABLE)
        total_rows += len(df_top)
        print(f"{len(df_top):,} rows ({time.time()-t0:.1f}s)")
    except Exception as e:
        print(f"ERROR: {e}")

cnt = spark.sql(f"SELECT COUNT(*) AS n FROM {TABLE}").collect()[0]["n"]
print(f"\n✅ {TABLE}: {cnt:,} rows")

# COMMAND ----------

# DBTITLE 1,Tahoe — tahoe_100m drug catalog (metadata sweep, no expression arrays)
from huggingface_hub import list_repo_files, hf_hub_download
import pandas as pd

TABLE = f"{PREFIX}.tahoe_100m"

# Build drug catalog: unique drug × cell_line × MOA rows.
# Sample every 17th shard (~200 of 3,388) — sufficient to capture all
# unique combinations without sweeping the full 241 GB dataset.
# Set SAMPLE_EVERY = 1 for exhaustive coverage.

META_COLS    = ["drug", "cell_line_id", "moa-fine", "canonical_smiles", "pubchem_cid", "plate"]
SAMPLE_EVERY = 17

all_files = sorted([
    f for f in list_repo_files(TAHOE_REPO, repo_type="dataset")
    if f.startswith("data/train-") and f.endswith(".parquet")
])
sampled = all_files[::SAMPLE_EVERY]
print(f"Total shards : {len(all_files):,}  |  Sampling: {len(sampled)} (every {SAMPLE_EVERY}th)")

catalog_dfs = []
for i, fname in enumerate(sampled):
    try:
        p = hf_hub_download(TAHOE_REPO, fname, repo_type="dataset")
        catalog_dfs.append(pd.read_parquet(p, columns=META_COLS))
    except Exception as e:
        print(f"  skip {fname}: {e}")
    if (i + 1) % 25 == 0:
        print(f"  {i+1}/{len(sampled)} shards ({sum(len(d) for d in catalog_dfs):,} rows)")

catalog = pd.concat(catalog_dfs, ignore_index=True)
dedup = (
    catalog
    .drop_duplicates(subset=["drug", "cell_line_id", "moa-fine"])
    .rename(columns={"moa-fine": "moa_fine"})
    .reset_index(drop=True)
)
print(f"\nUnique drug×cell_line rows : {len(dedup):,}")
print(f"Unique drugs               : {dedup['drug'].nunique():,}")
print(f"Unique cell lines           : {dedup['cell_line_id'].nunique():,}")
print(dedup.head(3).to_string())

spark.createDataFrame(dedup).write \
    .mode("overwrite").option("overwriteSchema", "true").saveAsTable(TABLE)

cnt      = spark.sql(f"SELECT COUNT(*) AS n FROM {TABLE}").collect()[0]["n"]
drug_cnt = spark.sql(f"SELECT COUNT(DISTINCT drug) AS n FROM {TABLE}").collect()[0]["n"]
print(f"\n✅ {TABLE}: {cnt:,} rows | {drug_cnt:,} unique drugs")

# COMMAND ----------

# DBTITLE 1,Grant SELECT on all new tables to app + SPs
SPS = [
    "2a8be974-109e-48a9-8d86-06e567ac7c3a",  # app SP
    "7b7e2188-3c15-4302-bf7c-770ef2e54c69",
    "1b52826b-7a90-488f-a605-c2468aaf3bf4",
]

new_tables = [
    f"{PREFIX}.crispr_atlas_gene_metadata",
    f"{PREFIX}.crispr_atlas_cell_metadata",
    f"{PREFIX}.crispr_atlas_ndd_genes",
    f"{PREFIX}.crispr_atlas_diff_expr",
    f"{PREFIX}.wholebrain_crispr_atlas",
    f"{PREFIX}.tahoe_100m_gene_vocab",
    f"{PREFIX}.tahoe_100m_clustered",
    f"{PREFIX}.tahoe_100m",
]

granted = errors = 0
for tbl in new_tables:
    for sp in SPS:
        try:
            spark.sql(f"GRANT SELECT ON TABLE {tbl} TO `{sp}`")
            granted += 1
        except Exception as e:
            if "not found" in str(e).lower():
                print(f"  ⚠️  {tbl.split('.')[-1]} not yet created — skip")
            else:
                print(f"  ❌ {tbl.split('.')[-1]} → {sp[:8]}: {e}")
                errors += 1

print(f"\n✅ Granted {granted} permissions ({errors} errors)")

# COMMAND ----------

# DBTITLE 1,Verify — row counts for all new tables
check_tables = [
    "crispr_atlas_gene_metadata",
    "crispr_atlas_cell_metadata",
    "crispr_atlas_ndd_genes",
    "crispr_atlas_diff_expr",
    "wholebrain_crispr_atlas",
    "tahoe_100m_gene_vocab",
    "tahoe_100m_clustered",
    "tahoe_100m",
]

print(f"{'Table':<38} {'Rows':>12}")
print("-" * 52)
for t in check_tables:
    fqn = f"{PREFIX}.{t}"
    try:
        cnt = spark.sql(f"SELECT COUNT(*) AS n FROM {fqn}").collect()[0]["n"]
        flag = "✅" if cnt > 0 else "⬜"
        print(f"{flag} {t:<36} {cnt:>12,}")
    except Exception as e:
        print(f"❌ {t:<36} NOT FOUND: {e}")

# COMMAND ----------

# DBTITLE 1,CRISPR Brain — complete load from CRISPRBrain.org
# CRISPRBrain.org human iPSC-neuron screens — complete load
# Current state: 23,040 rows / 2 screens.  This cell discovers and loads all available screens.
import os, sys, requests, io, time
import pandas as pd

os.environ.setdefault("NEUROPLEX_ENV", "prod")
for mod in list(sys.modules):
    if "ingestion" in mod or "neuroplex_config" in mod:
        del sys.modules[mod]
sys.path.insert(0, "/Workspace/Users/andrew_forman@eisai.com/neuro-crispr-mcp")
from ingestion.source_registry import TARGET_CATALOG, TARGET_SCHEMA
PREFIX_CB = f"{TARGET_CATALOG}.{TARGET_SCHEMA}"

TABLE = f"{PREFIX_CB}.crisprbrain_screens"
HDRS  = {"User-Agent": "NeuroPlex/1.0 (Eisai Discovery)", "Accept": "application/json"}
BASE  = "https://crisprbrain.org"

# ── 1. Current state ─────────────────────────────────────────────
existing_cnt = spark.sql(f"SELECT COUNT(*) AS n FROM {TABLE}").collect()[0]["n"]
existing_screens = [
    r[0] for r in spark.sql(
        f"SELECT DISTINCT screen_name FROM {TABLE} WHERE screen_name IS NOT NULL"
    ).collect()
]
print(f"Current: {existing_cnt:,} rows | {len(existing_screens)} screen(s): {existing_screens}")

# ── 2. Discover all available screens via API ─────────────────────
all_screens = []
for endpoint in ["/api/screens/", "/api/v1/screens/", "/api/"]:
    try:
        r = requests.get(f"{BASE}{endpoint}", headers=HDRS, timeout=20)
        if r.status_code == 200 and r.headers.get("Content-Type", "").startswith("application/json"):
            data = r.json()
            if isinstance(data, list):
                all_screens = data
            elif isinstance(data, dict):
                all_screens = data.get("results", data.get("screens", []))
            print(f"\nAPI {endpoint}: {len(all_screens)} screen(s) found")
            for s in all_screens[:5]:
                print(f"  {s}")
            break
        else:
            print(f"  {endpoint}: {r.status_code}")
    except Exception as e:
        print(f"  {endpoint}: {e}")

# ── 3. Bulk TSV / Excel fallback ──────────────────────────────────
new_rows = 0
bulk_loaded = False

for url in [
    f"{BASE}/api/bulk_download/?format=tsv",
    f"{BASE}/api/bulk/?format=tsv",
    f"{BASE}/download/bulk/all_phenotypes.tsv",
    f"{BASE}/static/downloads/all_phenotypes.tsv",
    f"{BASE}/download/all_phenotypes.tsv.gz",
]:
    try:
        print(f"\nTrying bulk URL: {url}")
        r = requests.get(url, headers={**HDRS, "Accept": "*/*"}, timeout=120, stream=True)
        if r.status_code != 200:
            print(f"  → {r.status_code}")
            continue
        ctype = r.headers.get("Content-Type", "")
        print(f"  → 200  Content-Type: {ctype}  Size: {r.headers.get('Content-Length','?')} bytes")
        raw = b"".join(r.iter_content(chunk_size=65536))
        sep = "\t" if (b"\t" in raw[:2000] or "tsv" in ctype) else ","
        df_bulk = pd.read_csv(io.BytesIO(raw), sep=sep, low_memory=False)
        print(f"  Parsed: {df_bulk.shape}  cols: {list(df_bulk.columns[:8])}")
        # Filter to rows with phenotype data (skip guide-library sheets)
        if "has_phenotype" in df_bulk.columns:
            df_bulk = df_bulk[df_bulk["has_phenotype"].astype(str).str.lower().isin(["true", "1", "yes"])]
        # Append only screens not already in the table
        scr_col = next((c for c in df_bulk.columns if "screen" in c.lower()), None)
        if scr_col and existing_screens:
            df_new = df_bulk[~df_bulk[scr_col].isin(existing_screens)]
        else:
            df_new = df_bulk
        if len(df_new) > 0:
            sdf = spark.createDataFrame(df_new.astype(str))
            sdf.write.mode("append").option("mergeSchema", "true").saveAsTable(TABLE)
            new_rows += len(df_new)
            print(f"  ✅ Appended {len(df_new):,} new rows")
        bulk_loaded = True
        break
    except Exception as e:
        print(f"  {e}")

# ── 4. Per-screen Excel download (if bulk not available) ──────────
if not bulk_loaded and all_screens:
    print("\nBulk download unavailable — fetching screens individually…")
    %pip install -q openpyxl
    import openpyxl

    def screen_id(s):
        return s.get("id", s.get("name", s.get("slug", str(s)))) if isinstance(s, dict) else str(s)

    for screen in all_screens:
        sid = screen_id(screen)
        if sid in existing_screens:
            print(f"  {sid}: already loaded — skip")
            continue
        # Try known URL patterns for individual screen downloads
        for tmpl in [
            f"{BASE}/api/screens/{sid}/download/",
            f"{BASE}/static/downloads/{sid}.xlsx",
            f"{BASE}/download/{sid}.xlsx",
        ]:
            try:
                r = requests.get(tmpl, headers={**HDRS, "Accept": "*/*"}, timeout=60)
                if r.status_code != 200:
                    continue
                print(f"  {sid}: fetched from {tmpl} ({len(r.content):,} bytes)")
                wb = openpyxl.load_workbook(io.BytesIO(r.content), read_only=True, data_only=True)
                rows_this = 0
                for sheet in wb.sheetnames:
                    ws = wb[sheet]
                    data = list(ws.values)
                    if not data:
                        continue
                    headers = [str(h) for h in data[0]]
                    if "has_phenotype" not in headers and "phenotype" not in " ".join(headers).lower():
                        continue   # guide-library or non-phenotype sheet
                    df_sheet = pd.DataFrame(data[1:], columns=headers)
                    df_sheet["screen_name"] = sid
                    df_sheet["source_sheet"] = sheet
                    sdf = spark.createDataFrame(df_sheet.astype(str))
                    sdf.write.mode("append").option("mergeSchema", "true").saveAsTable(TABLE)
                    rows_this += len(df_sheet)
                new_rows += rows_this
                print(f"    → {rows_this:,} rows written")
                break
            except Exception as e:
                pass

# ── 5. Final summary ──────────────────────────────────────────────
final_cnt = spark.sql(f"SELECT COUNT(*) AS n FROM {TABLE}").collect()[0]["n"]
final_screens = [
    r[0] for r in spark.sql(
        f"SELECT DISTINCT screen_name FROM {TABLE} WHERE screen_name IS NOT NULL"
    ).collect()
]
print(f"\n{'─'*55}")
print(f"✅ crisprbrain_screens final: {final_cnt:,} rows")
print(f"   Screens loaded ({len(final_screens)}): {sorted(final_screens)}")
print(f"   New rows added: {final_cnt - existing_cnt:,}")

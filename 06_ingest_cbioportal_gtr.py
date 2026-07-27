# Databricks notebook source
# /// script
# [tool.databricks.environment]
# environment_version = "5"
# ///
# DBTITLE 1,Config — environment, gene panel, API endpoints
import os, sys, json, time, logging
import requests
from collections import defaultdict
from pyspark.sql import functions as F
from pyspark.sql.types import StructType, StructField, StringType

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("cbio_gtr_etl")

os.environ["NEUROPLEX_ENV"] = "val"
sys.path.insert(0, "/Workspace/Users/andrew_forman@eisai.com/neuro-crispr-mcp")
for mod in [k for k in sys.modules if k.startswith("config")]:
    del sys.modules[mod]
from config.neuroplex_config import load_config

CFG = load_config()
CAT, SCH = CFG.catalog, CFG.schema

TABLE_CBIO = f"{CAT}.{SCH}.neuroplex_cbioportal"
TABLE_GTR  = f"{CAT}.{SCH}.neuroplex_ncbi_gtr"

GENE_PANEL = [
    "PSEN1", "PSEN2", "APP",    "APOE",   "TREM2",  "BIN1",
    "CLU",   "CR1",   "PICALM", "ABCA7",  "SORL1",
    "SNCA",  "LRRK2", "PINK1",  "PRKN",   "PARK7",  "GBA1",  "VPS35",
    "SOD1",  "TARDBP","FUS",    "C9orf72", "TBK1",  "GRN",
    "MAPT",  "HTT",
    "HCRT",  "HCRTR1","HCRTR2",
    "MECP2", "FMR1",
]

CBIO_API   = "https://www.cbioportal.org/api"
EUTILS_API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# NCBI recommends <=3 req/s without an API key
# Set NCBI_API_KEY env var to allow 10 req/s
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")
CBIO_DELAY  = 0.5   # cBioPortal is generally permissive
NCBI_DELAY  = 0.4   # ~2.5 req/s — safe without an API key

print(f"✅ Environment : {CFG.environment}")
print(f"   Catalog     : {CAT}.{SCH}")
print(f"   cBioPortal  : {TABLE_CBIO}")
print(f"   NCBI GTR    : {TABLE_GTR}")
print(f"   Gene panel  : {len(GENE_PANEL)} genes")
if NCBI_API_KEY:
    print(f"   NCBI key    : ...{NCBI_API_KEY[-4:]} (10 req/s allowed)")

# COMMAND ----------

# DBTITLE 1,cBioPortal REST API helpers
# cBioPortal REST API — somatic mutation landscape across TCGA PanCancer Atlas 2018
# Pipeline:
#   1. GET /genes/{symbol} → entrezGeneId
#   2. POST /mutations/fetch → somatic mutations in TCGA PanCan Atlas
#   3. Aggregate by proteinChange + mutationType → top N rows
#
# neuroplex_cbioportal payload shape matches neuroplex_query.py:
#   payload:mutation_type, payload:protein_change, payload:cancer_study

# TCGA PanCancer Atlas 2018 — most comprehensive cross-cancer somatic mutation dataset
# PanCan Atlas is split into 33 per-cancer studies; aggregate across 8 representative ones
TCGA_STUDIES = [
    ("gbm_tcga_pan_can_atlas_2018",      "Glioblastoma (TCGA)"),
    ("lgg_tcga_pan_can_atlas_2018",      "Lower Grade Glioma (TCGA)"),
    ("laml_tcga_pan_can_atlas_2018",     "Acute Myeloid Leukemia (TCGA)"),
    ("brca_tcga_pan_can_atlas_2018",     "Breast Carcinoma (TCGA)"),
    ("luad_tcga_pan_can_atlas_2018",     "Lung Adenocarcinoma (TCGA)"),
    ("coadread_tcga_pan_can_atlas_2018", "Colorectal Adenocarcinoma (TCGA)"),
    ("skcm_tcga_pan_can_atlas_2018",     "Skin Melanoma (TCGA)"),
    ("stad_tcga_pan_can_atlas_2018",     "Stomach Adenocarcinoma (TCGA)"),
]


def _cbio_get(path: str, params: dict | None = None, retries: int = 3) -> list | dict | None:
    url = f"{CBIO_API}/{path.lstrip('/')}"
    headers = {"Accept": "application/json"}
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, headers=headers, timeout=30)
            if r.status_code == 200:
                return r.json()
            elif r.status_code == 429:
                wait = 2 ** (attempt + 1)
                log.warning("cBioPortal rate limited, waiting %ss", wait)
                time.sleep(wait)
            else:
                log.warning("cBioPortal HTTP %s: %s", r.status_code, r.text[:200])
                return None
        except Exception as exc:
            log.warning("cBioPortal error (attempt %d): %s", attempt + 1, exc)
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
    return None


def _cbio_post(path: str, body: dict, retries: int = 3) -> list | dict | None:
    url = f"{CBIO_API}/{path.lstrip('/')}"
    headers = {"Accept": "application/json", "Content-Type": "application/json"}
    for attempt in range(retries):
        try:
            r = requests.post(url, json=body, headers=headers, timeout=60)
            if r.status_code == 200:
                return r.json()
            elif r.status_code == 429:
                wait = 2 ** (attempt + 1)
                time.sleep(wait)
            else:
                log.warning("cBioPortal POST HTTP %s: %s", r.status_code, r.text[:200])
                return None
        except Exception as exc:
            log.warning("cBioPortal POST error (attempt %d): %s", attempt + 1, exc)
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
    return None


def resolve_cbio_gene(gene_symbol: str) -> int | None:
    """Resolve gene symbol → Entrez Gene ID via cBioPortal."""
    result = _cbio_get(f"genes/{gene_symbol}", params={"geneIdType": "HUGO_GENE_SYMBOL"})
    if isinstance(result, dict):
        return result.get("entrezGeneId")
    return None


def fetch_cbio_mutations(gene_symbol: str, max_rows: int = 50) -> list[dict]:
    """Fetch somatic mutations across 8 representative TCGA PanCan Atlas 2018 studies."""
    entrez_id = resolve_cbio_gene(gene_symbol)
    if not entrez_id:
        log.warning("Could not resolve Entrez ID for %s", gene_symbol)
        return []

    agg: dict[tuple, dict] = {}
    for study_id, study_label in TCGA_STUDIES:
        profile_id  = f"{study_id}_mutations"
        sample_list = f"{study_id}_all"
        mutations = _cbio_post(
            f"molecular-profiles/{profile_id}/mutations/fetch?projection=SUMMARY",
            body={"sampleListId": sample_list, "entrezGeneIds": [entrez_id]},
        )
        if not mutations:
            time.sleep(CBIO_DELAY * 0.5)
            continue
        for m in (mutations if isinstance(mutations, list) else []):
            key = (m.get("mutationType") or "Unknown", m.get("proteinChange") or "")
            if key not in agg:
                agg[key] = {"count": 0, "cancer_types": set()}
            agg[key]["count"] += 1
            agg[key]["cancer_types"].add(study_label)
        time.sleep(CBIO_DELAY * 0.5)

    sorted_muts = sorted(agg.items(), key=lambda x: -x[1]["count"])
    records = []
    for (mut_type, protein_chg), info in sorted_muts[:max_rows]:
        cancer_str = "; ".join(sorted(info["cancer_types"])[:3])
        records.append({
            "gene_symbol": gene_symbol,
            "title":   protein_chg or mut_type,
            "disease": cancer_str,
            "summary": f"{mut_type} ({info['count']} samples, TCGA Pan-Cancer)",
            "payload": {
                "mutation_type":  mut_type,
                "protein_change": protein_chg,
                "cancer_study":   "TCGA PanCancer Atlas 2018",
                "sample_count":   info["count"],
                "cancer_types":   sorted(info["cancer_types"]),
                "entrez_id":      entrez_id,
            },
        })
    return records


print("✅ cBioPortal helpers defined")

# COMMAND ----------

# DBTITLE 1,NCBI GTR E-utilities helpers
# NCBI GTR via E-utilities (public, API key optional but recommended for >3 req/s)
# Pipeline per gene:
#   1. esearch db=gene term={symbol}[gene_name]+AND+Homo+sapiens[organism]
#      → NCBI Gene ID (entrez_gene_id)
#   2. esearch db=gtr term={ncbi_gene_id}[geneid]
#      → list of GTR test UIDs
#   3. esummary db=gtr id={uids}
#      → test name, condition, test type, lab name
#
# neuroplex_ncbi_gtr payload shape matches neuroplex_query.py:
#   payload:test_type, payload:lab_name

def _ncbi_get(endpoint: str, params: dict, retries: int = 3) -> dict | None:
    params.setdefault("retmode", "json")
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    url = f"{EUTILS_API}/{endpoint}"
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=30,
                             headers={"Accept": "application/json"})
            if r.status_code == 200:
                return r.json()
            elif r.status_code == 429:
                wait = 2 ** (attempt + 1)
                log.warning("NCBI rate limited, waiting %ss", wait)
                time.sleep(wait)
            else:
                log.warning("NCBI HTTP %s: %s", r.status_code, r.text[:200])
                return None
        except Exception as exc:
            log.warning("NCBI error (attempt %d): %s", attempt + 1, exc)
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
    return None


def _ncbi_gene_id(gene_symbol: str) -> str | None:
    """Resolve gene symbol → NCBI Gene ID (Entrez) via Gene database."""
    # Use [gene] field tag (not [gene_name]) and let requests encode spaces properly.
    # The [gene] field searches the official gene symbol / name in NCBI Gene.
    result = _ncbi_get("esearch.fcgi", {
        "db":     "gene",
        "term":   f"{gene_symbol}[gene] AND Homo sapiens[organism] AND alive[property]",
        "retmax": "3",
        "sort":   "relevance",
    })
    ids = ((result or {}).get("esearchresult") or {}).get("idlist") or []
    return ids[0] if ids else None


def fetch_gtr_tests(gene_symbol: str, max_tests: int = 20) -> list[dict]:
    """Fetch NCBI GTR genetic test listings for a gene."""
    ncbi_id = _ncbi_gene_id(gene_symbol)
    if not ncbi_id:
        log.warning("No NCBI Gene ID for %s", gene_symbol)
        return []
    time.sleep(NCBI_DELAY)

    # Search GTR tests linked to this gene
    search = _ncbi_get("esearch.fcgi", {
        "db": "gtr",
        "term": f"{ncbi_id}[geneid]",
        "retmax": str(max_tests),
    })
    uids = ((search or {}).get("esearchresult") or {}).get("idlist") or []
    if not uids:
        return []
    time.sleep(NCBI_DELAY)

    # Fetch test summaries
    summary = _ncbi_get("esummary.fcgi", {
        "db": "gtr",
        "id": ",".join(uids),
    })
    result_block = (summary or {}).get("result") or {}
    if not result_block:
        return []

    records = []
    for uid in uids:
        test = result_block.get(uid) or result_block.get(str(uid))
        if not test or uid == "uids":
            continue

        test_name    = test.get("testname") or test.get("title") or f"GTR{uid}"
        test_type    = test.get("testtype") or ""
        lab_name     = test.get("labname") or ""
        gtr_id       = test.get("gtid") or uid

        # Conditions tested: can be a list or a comma-joined string
        conditions   = test.get("conditionname") or ""
        if isinstance(conditions, list):
            conditions = "; ".join(c.get("name", "") if isinstance(c, dict) else str(c)
                                   for c in conditions if c)[:300]
        else:
            conditions = str(conditions)[:300]

        records.append({
            "gene_symbol": gene_symbol,
            "title": test_name[:200],
            "summary": f"{test_type} test for {conditions[:100]}" if conditions else test_name[:200],
            "disease": conditions[:300],
            "payload": {
                "test_type":     test_type,
                "lab_name":      lab_name,
                "gtr_id":        gtr_id,
                "ncbi_gene_id":  ncbi_id,
                "conditions":    conditions,
                "test_name":     test_name,
            },
        })
    return records


print("✅ NCBI GTR helpers defined")

# COMMAND ----------

# DBTITLE 1,Connectivity probe — PSEN1 + SOD1
# cBioPortal probe: PSEN1 (germline AD gene → few somatic mutations)
#                   SOD1  (ALS gene → expected somatic mutations in some cancers)
for gene in ["SOD1", "PSEN1"]:
    eid = resolve_cbio_gene(gene)
    print(f"{gene}: Entrez ID = {eid}")
    muts = fetch_cbio_mutations(gene, max_rows=5)
    if muts:
        for m in muts:
            print(f"  {m['payload']['mutation_type']:30s}  {m['payload']['protein_change']:15s}  n={m['payload']['sample_count']}")
    else:
        print(f"  (no somatic mutations in TCGA Pan-Can Atlas — expected for germline genes)")
    print()
    time.sleep(CBIO_DELAY)

print("\n--- GTR probe: PSEN1 ---")
gtr = fetch_gtr_tests("PSEN1", max_tests=3)
for t in gtr:
    print(f"  {t['title'][:60]:60s}  [{t['payload']['test_type']}]  {t['payload']['lab_name'][:30]}")
if not gtr:
    print("  (no tests found)")

print("\n✅ Connectivity probes complete")

# COMMAND ----------

# DBTITLE 1,Fetch cBioPortal mutations for all panel genes
cbio_records = []
cbio_errors  = []

for i, gene in enumerate(GENE_PANEL):
    print(f"[{i+1:02d}/{len(GENE_PANEL)}] {gene} ...", end=" ")
    try:
        recs = fetch_cbio_mutations(gene, max_rows=50)
        cbio_records.extend(recs)
        if recs:
            total = sum(r["payload"]["sample_count"] for r in recs)
            print(f"✅  {len(recs)} unique mutations, {total} total samples")
        else:
            print("⚠️  no somatic mutations in TCGA Pan-Can Atlas")
    except Exception as exc:
        cbio_errors.append((gene, str(exc)))
        print(f"❌  {exc}")
    time.sleep(CBIO_DELAY)

print(f"\n———")
print(f"cBioPortal records : {len(cbio_records):,}")
print(f"Genes with mutations: {len(set(r['gene_symbol'] for r in cbio_records))}/{len(GENE_PANEL)}")
if cbio_errors:
    print(f"Errors: {cbio_errors}")

# COMMAND ----------

# DBTITLE 1,Fetch NCBI GTR tests for all panel genes
gtr_records = []
gtr_errors  = []

for i, gene in enumerate(GENE_PANEL):
    print(f"[{i+1:02d}/{len(GENE_PANEL)}] {gene} ...", end=" ")
    try:
        recs = fetch_gtr_tests(gene, max_tests=20)
        gtr_records.extend(recs)
        if recs:
            diseases = "; ".join({r['disease'][:30] for r in recs if r['disease']})[:80]
            print(f"✅  {len(recs)} tests  ({diseases}...)")
        else:
            print("⚠️  no GTR tests found")
    except Exception as exc:
        gtr_errors.append((gene, str(exc)))
        print(f"❌  {exc}")
    time.sleep(NCBI_DELAY)

print(f"\n———")
print(f"GTR records   : {len(gtr_records):,}")
print(f"Genes covered : {len(set(r['gene_symbol'] for r in gtr_records))}/{len(GENE_PANEL)}")
if gtr_errors:
    print(f"Errors: {gtr_errors}")

# COMMAND ----------

# DBTITLE 1,Write neuroplex_cbioportal + neuroplex_ncbi_gtr to Delta
def _write_records(records: list[dict], table: str, str_cols: list[str]):
    """Write records with string + VARIANT payload columns to a Delta table."""
    if not records:
        print(f"  ⚠️  No records to write to {table} — writing empty placeholder")
        # Keep the table schema intact but empty
        spark.sql(f"DELETE FROM {table} WHERE 1=1")
        return

    rows = [tuple(r.get(c, "") for c in str_cols) + (json.dumps(r.get("payload", {}), default=str),)
            for r in records]

    schema = StructType(
        [StructField(c, StringType(), True) for c in str_cols]
        + [StructField("payload_str", StringType(), True)]
    )
    df_raw = spark.createDataFrame(rows, schema)
    df_final = df_raw.select(
        *[F.col(c) for c in str_cols],
        F.expr("parse_json(payload_str)").alias("payload"),
    )
    n = df_final.count()
    print(f"  Rows to write: {n:,}")
    (
        df_final.write
        .format("delta")
        .mode("overwrite")
        .option("overwriteSchema", "false")
        .saveAsTable(table)
    )
    print(f"  ✅ Written to {table}")


print("Writing neuroplex_cbioportal ...")
_write_records(
    cbio_records, TABLE_CBIO,
    str_cols=["gene_symbol", "title", "disease", "summary"]
)

print("\nWriting neuroplex_ncbi_gtr ...")
_write_records(
    gtr_records, TABLE_GTR,
    str_cols=["gene_symbol", "title", "summary", "disease"]
)

# COMMAND ----------

# DBTITLE 1,Validation — cBioPortal spot-check
# MAGIC %sql
# MAGIC -- Mutation landscape by gene: expect ALS/cancer genes (SOD1, TARDBP, FUS) to have high counts
# MAGIC SELECT
# MAGIC     gene_symbol,
# MAGIC     COUNT(*)                                        AS unique_mutations,
# MAGIC     SUM(payload:sample_count::INT)                  AS total_samples_mutated,
# MAGIC     FIRST(payload:mutation_type::STRING)            AS top_mutation_type,
# MAGIC     FIRST(payload:protein_change::STRING)           AS top_protein_change
# MAGIC FROM dhbl_discovery_us_val.genesis_schema.neuroplex_cbioportal
# MAGIC GROUP BY gene_symbol
# MAGIC ORDER BY total_samples_mutated DESC

# COMMAND ----------

# DBTITLE 1,Validation — GTR spot-check
# MAGIC %sql
# MAGIC -- PSEN1 GTR tests: expect Alzheimer panel tests, predictive testing, diagnostic tests
# MAGIC SELECT
# MAGIC     gene_symbol,
# MAGIC     title,
# MAGIC     payload:test_type::STRING   AS test_type,
# MAGIC     payload:lab_name::STRING    AS lab_name,
# MAGIC     LEFT(disease, 100)          AS condition
# MAGIC FROM dhbl_discovery_us_val.genesis_schema.neuroplex_ncbi_gtr
# MAGIC WHERE UPPER(gene_symbol) = 'PSEN1'
# MAGIC ORDER BY title
# MAGIC LIMIT 15

# COMMAND ----------

# DBTITLE 1,GBA1 patch — append to existing tables
# Targeted GBA1 append — runs independently of cells 5-6.
# Helper functions must be in memory (run cells 2-3 if kernel was restarted).

def _append_gba1(records, table, str_cols):
    if not records:
        print(f"  ⚠️  No records returned for GBA1 → {table}")
        return
    rows = [
        tuple(r.get(c, "") for c in str_cols)
        + (json.dumps(r.get("payload", {}), default=str),)
        for r in records
    ]
    schema = StructType(
        [StructField(c, StringType(), True) for c in str_cols]
        + [StructField("payload_str", StringType(), True)]
    )
    df_raw = spark.createDataFrame(rows, schema)
    df_final = df_raw.select(
        *[F.col(c) for c in str_cols],
        F.expr("parse_json(payload_str)").alias("payload"),
    )
    n = df_final.count()
    df_final.write.format("delta").mode("append").saveAsTable(table)
    print(f"  ✅ Appended {n} GBA1 rows → {table}")


print("Fetching GBA1 → cBioPortal ...")
gba1_cbio = fetch_cbio_mutations("GBA1", max_rows=50)
if gba1_cbio:
    total = sum(r["payload"]["sample_count"] for r in gba1_cbio)
    print(f"  {len(gba1_cbio)} unique mutations, {total} total samples")
_append_gba1(gba1_cbio, TABLE_CBIO, ["gene_symbol", "title", "disease", "summary"])

time.sleep(CBIO_DELAY)

print("\nFetching GBA1 → NCBI GTR ...")
gba1_gtr = fetch_gtr_tests("GBA1", max_tests=20)
if gba1_gtr:
    print(f"  {len(gba1_gtr)} tests found")
_append_gba1(gba1_gtr, TABLE_GTR, ["gene_symbol", "title", "summary", "disease"])

# COMMAND ----------

# DBTITLE 1,Deploy neuro-crispr-mcp app
from databricks.sdk import WorkspaceClient
from databricks.sdk.service.apps import AppDeployment

w = WorkspaceClient()
deployment = w.apps.deploy(
    app_name="neuro-crispr-mcp",
    app_deployment=AppDeployment(
        source_code_path="/Workspace/Users/andrew_forman@eisai.com/neuro-crispr-mcp"
    ),
).result()  # .result() waits for deploy to complete

print(f"Deployment ID : {deployment.deployment_id}")
print(f"Source path   : {deployment.source_code_path}")
print(f"Status        : {deployment.status.state if deployment.status else 'unknown'}")
print(f"Status msg    : {deployment.status.message if deployment.status else ''}")
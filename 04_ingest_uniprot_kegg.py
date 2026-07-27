# Databricks notebook source
# /// script
# [tool.databricks.environment]
# environment_version = "5"
# ///
# DBTITLE 1,Config — environment, gene panel, API endpoints
import os, sys, json, time, logging
import requests
from pyspark.sql import functions as F
from pyspark.sql.types import StructType, StructField, StringType

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("uniprot_kegg_etl")

os.environ["NEUROPLEX_ENV"] = "val"
sys.path.insert(0, "/Workspace/Users/andrew_forman@eisai.com/neuro-crispr-mcp")
for mod in [k for k in sys.modules if k.startswith("config")]:
    del sys.modules[mod]
from config.neuroplex_config import load_config

CFG = load_config()
CAT = CFG.catalog
SCH = CFG.schema

TABLE_UNIPROT = f"{CAT}.{SCH}.neuroplex_uniprot"
TABLE_KEGG    = f"{CAT}.{SCH}.neuroplex_kegg"

print(f"✅ Environment : {CFG.environment}")
print(f"   Catalog     : {CAT}.{SCH}")
print(f"   Targets     : {TABLE_UNIPROT}")
print(f"                 {TABLE_KEGG}")

# ── NeuroPlex gene panel (same 31 genes as gnomAD ETL) ──────────────────────
GENE_PANEL = [
    "PSEN1", "PSEN2", "APP",   "APOE",  "TREM2",  "BIN1",
    "CLU",   "CR1",   "PICALM","ABCA7", "SORL1",
    "SNCA",  "LRRK2", "PINK1", "PRKN",  "PARK7",  "GBA1",  "VPS35",
    "SOD1",  "TARDBP","FUS",   "C9orf72","TBK1",  "GRN",
    "MAPT",  "HTT",
    "HCRT",  "HCRTR1","HCRTR2",
    "MECP2", "FMR1",
]

# API base URLs (both public, no auth required)
UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"
KEGG_API    = "https://rest.kegg.jp"
REQUEST_DELAY = 1.0  # seconds between calls

print(f"   Gene panel  : {len(GENE_PANEL)} genes")

# COMMAND ----------

# DBTITLE 1,UniProt REST API helpers
# UniProt REST API v2 — reviewed (SwissProt) human entries only.
# Fields returned: accession, gene names, protein name, function, subcellular location.
# One primary entry per gene (canonical isoform); append isoforms if needed later.

UNIPROT_FIELDS = ",".join([
    "accession",
    "gene_names",
    "protein_name",
    "cc_function",
    "cc_subcellular_location",
    "keyword",
    "organism_name",
])


def _http_get(url: str, params: dict | None = None, retries: int = 3) -> requests.Response | None:
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=30,
                             headers={"Accept": "application/json"})
            if r.status_code == 200:
                return r
            elif r.status_code == 429:
                wait = 2 ** (attempt + 1)
                log.warning("Rate limited, waiting %ss", wait)
                time.sleep(wait)
            else:
                log.warning("HTTP %s: %s", r.status_code, r.text[:200])
                return None
        except Exception as exc:
            log.warning("Request error (attempt %d): %s", attempt + 1, exc)
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
    return None


def _extract_comment(comments: list, comment_type: str) -> str:
    """Pull free-text from a UniProt comments block of a given type."""
    for c in (comments or []):
        if c.get("commentType") == comment_type:
            texts = c.get("texts") or []
            if texts:
                return texts[0].get("value", "")
            # subcellular location uses a nested structure
            locs = c.get("subcellularLocations") or []
            if locs:
                parts = []
                for loc in locs:
                    loc_name = (loc.get("location") or {}).get("value")
                    if loc_name:
                        parts.append(loc_name)
                return "; ".join(parts)
    return ""


def fetch_uniprot(gene_symbol: str) -> dict | None:
    """Fetch the canonical SwissProt entry for a human gene. Returns a record dict or None."""
    r = _http_get(UNIPROT_API, params={
        "query":  f'gene_exact:"{gene_symbol}" AND organism_id:9606 AND reviewed:true',
        "fields": UNIPROT_FIELDS,
        "format": "json",
        "size":   "1",
    })
    if not r:
        return None

    results = r.json().get("results") or []
    if not results:
        return None

    entry = results[0]
    acc   = entry.get("primaryAccession", "")
    comments = entry.get("comments") or []

    # Protein full name
    prot_desc = entry.get("proteinDescription") or {}
    rec_name  = prot_desc.get("recommendedName") or {}
    prot_name = (rec_name.get("fullName") or {}).get("value", "")

    function_text = _extract_comment(comments, "FUNCTION")
    location_text = _extract_comment(comments, "SUBCELLULAR LOCATION")

    # Keywords (comma-joined)
    keywords = ", ".join(
        (kw.get("name") or kw.get("value") or "")
        for kw in (entry.get("keywords") or [])
    )[:500]  # truncate for storage

    return {
        "gene_symbol": gene_symbol,
        "title": f"UniProt: {gene_symbol} ({acc})",
        "summary": function_text[:500] if function_text else prot_name,
        "payload": {
            "uniprot_acc": acc,
            "protein_name": prot_name,
            "function": function_text,
            "subcellular_location": location_text,
            "keywords": keywords,
        },
    }


print("✅ UniProt helpers defined")

# COMMAND ----------

# DBTITLE 1,KEGG REST API helpers
# KEGG REST API (text/plain responses, tab-delimited)
# Pipeline per gene:
#   1. find/genes/{gene_symbol} → filter for hsa: to get KEGG human gene IDs
#   2. link/pathway/{kegg_gene_id} → get pathway IDs (hsa05010, etc.)
#   3. Batch list/{ids} → get human-readable pathway names


def _kegg_get(path: str, retries: int = 3) -> str | None:
    """GET a KEGG REST endpoint, return response text or None."""
    url = f"{KEGG_API}/{path.lstrip('/')}"
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=30, headers={"Accept": "text/plain"})
            if r.status_code == 200:
                return r.text
            elif r.status_code == 429:
                wait = 2 ** (attempt + 1)
                log.warning("KEGG rate limited, waiting %ss", wait)
                time.sleep(wait)
            else:
                log.warning("KEGG HTTP %s for %s", r.status_code, url)
                return None
        except Exception as exc:
            log.warning("KEGG request error (attempt %d): %s", attempt + 1, exc)
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
    return None


def _kegg_gene_ids(gene_symbol: str) -> list[str]:
    """Return KEGG human gene IDs (hsa:XXXXX) for a gene symbol."""
    text = _kegg_get(f"find/genes/{gene_symbol}")
    if not text:
        return []
    ids = []
    for line in text.strip().splitlines():
        parts = line.split("\t")
        if not parts:
            continue
        kid = parts[0].strip()   # e.g. "hsa:5663"
        if not kid.startswith("hsa:"):
            continue
        # Confirm gene symbol matches (KEGG find is substring, not exact)
        desc = parts[1].upper() if len(parts) > 1 else ""
        syms = desc.split(" ")[0].split(";")[0]  # first token before space/semicolon
        if gene_symbol.upper() in syms.split(","):
            ids.append(kid)
    return ids


def _kegg_pathway_ids(kegg_gene_id: str) -> list[str]:
    """Return pathway IDs for a KEGG gene ID."""
    text = _kegg_get(f"link/pathway/{kegg_gene_id}")
    if not text:
        return []
    paths = []
    for line in text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            pid = parts[1].strip()   # e.g. "path:hsa05010"
            # Keep only human pathways (hsa)
            if ":hsa" in pid:
                paths.append(pid.replace("path:", ""))   # -> "hsa05010"
    return paths


# Module-level cache: fetched once, reused for all genes
_KEGG_HUMAN_PATHWAYS: dict[str, str] = {}


def _load_kegg_human_pathways() -> dict[str, str]:
    """Fetch all human KEGG pathways in one call and cache them.
    Returns dict mapping pathway_id (e.g. 'hsa05010') → clean name."""
    global _KEGG_HUMAN_PATHWAYS
    if _KEGG_HUMAN_PATHWAYS:
        return _KEGG_HUMAN_PATHWAYS
    text = _kegg_get("list/pathway/hsa")
    if not text:
        log.warning("Could not fetch KEGG human pathway list")
        return {}
    for line in text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            pid   = parts[0].replace("path:", "").strip()   # hsa05010
            pname = parts[1].replace(" - Homo sapiens (human)", "").strip()
            _KEGG_HUMAN_PATHWAYS[pid] = pname
    print(f"  Loaded {len(_KEGG_HUMAN_PATHWAYS)} human KEGG pathways")
    return _KEGG_HUMAN_PATHWAYS


def _kegg_pathway_names(pathway_ids: list[str]) -> dict[str, str]:
    """Look up pathway names from the cached human pathway map."""
    all_pathways = _load_kegg_human_pathways()
    return {pid: all_pathways.get(pid, pid) for pid in pathway_ids}


def fetch_kegg_pathways(gene_symbol: str) -> list[dict]:
    """Fetch all human KEGG pathway memberships for a gene. Returns list of record dicts."""
    gene_ids = _kegg_gene_ids(gene_symbol)
    if not gene_ids:
        return []

    all_pathway_ids = []
    for gid in gene_ids:
        all_pathway_ids.extend(_kegg_pathway_ids(gid))
        time.sleep(0.3)

    # Deduplicate
    all_pathway_ids = list(dict.fromkeys(all_pathway_ids))
    if not all_pathway_ids:
        return []

    pathway_names = _kegg_pathway_names(all_pathway_ids)

    records = []
    for pid in all_pathway_ids:
        pname = pathway_names.get(pid, pid)   # fallback to ID if name missing
        # Strip " - Homo sapiens (human)" suffix for cleaner display
        pname_clean = pname.replace(" - Homo sapiens (human)", "").strip()
        records.append({
            "gene_symbol": gene_symbol,
            "title": pname_clean,
            "summary": f"KEGG pathway {pid}: {pname_clean}",
            "payload": {
                "pathway_id":   pid,
                "pathway_name": pname_clean,
                "kegg_gene_ids": gene_ids,
            },
        })
    return records


print("✅ KEGG helpers defined")

# COMMAND ----------

# DBTITLE 1,API connectivity test — PSEN1 single-gene probe
# Quick sanity check before running the full panel.
print("=== UniProt probe: PSEN1 ===")
u = fetch_uniprot("PSEN1")
if u:
    print(f"  Accession  : {u['payload']['uniprot_acc']}")
    print(f"  Protein    : {u['payload']['protein_name']}")
    print(f"  Function   : {u['payload']['function'][:120]}...")
    print(f"  Location   : {u['payload']['subcellular_location']}")
else:
    print("  ⚠️  No result — check API")

time.sleep(REQUEST_DELAY)

print("\n=== KEGG probe: PSEN1 ===")
kr = fetch_kegg_pathways("PSEN1")
if kr:
    for p in kr[:5]:
        print(f"  {p['payload']['pathway_id']:12s}  {p['title']}")
else:
    print("  ⚠️  No pathways returned — check API")

# COMMAND ----------

# DBTITLE 1,Fetch UniProt for all panel genes
uniprot_records = []
uniprot_errors  = []

for i, gene in enumerate(GENE_PANEL):
    print(f"[{i+1:02d}/{len(GENE_PANEL)}] {gene} ...", end=" ")
    try:
        rec = fetch_uniprot(gene)
        if rec:
            uniprot_records.append(rec)
            print(f"✅  {rec['payload']['uniprot_acc']}  {rec['payload']['protein_name'][:50]}")
        else:
            print("⚠️  not found")
    except Exception as exc:
        uniprot_errors.append((gene, str(exc)))
        print(f"❌  {exc}")
    time.sleep(REQUEST_DELAY)

print(f"\n———")
print(f"UniProt records  : {len(uniprot_records)}/{len(GENE_PANEL)}")
if uniprot_errors:
    print(f"Errors ({len(uniprot_errors)}): {uniprot_errors}")

# COMMAND ----------

# DBTITLE 1,Fetch KEGG pathways for all panel genes
kegg_records = []
kegg_errors  = []

for i, gene in enumerate(GENE_PANEL):
    print(f"[{i+1:02d}/{len(GENE_PANEL)}] {gene} ...", end=" ")
    try:
        recs = fetch_kegg_pathways(gene)
        kegg_records.extend(recs)
        if recs:
            paths = ", ".join(r["payload"]["pathway_id"] for r in recs[:3])
            print(f"✅  {len(recs)} pathways ({paths}{' ...' if len(recs) > 3 else ''})")
        else:
            print("⚠️  no pathways found")
    except Exception as exc:
        kegg_errors.append((gene, str(exc)))
        print(f"❌  {exc}")
    time.sleep(REQUEST_DELAY)

print(f"\n———")
print(f"KEGG pathway records: {len(kegg_records)}  across {len(GENE_PANEL)} genes")
unique_paths = len(set(r['payload']['pathway_id'] for r in kegg_records))
print(f"Unique pathways      : {unique_paths}")
if kegg_errors:
    print(f"Errors ({len(kegg_errors)}): {kegg_errors}")

# COMMAND ----------

# DBTITLE 1,Write neuroplex_uniprot to Delta
def _write_records(records: list[dict], table: str, extra_cols: list[str] | None = None):
    """Write gene/title/summary/payload records to a Delta table using parse_json()."""
    if not records:
        raise RuntimeError(f"No records to write to {table}")

    base_cols = ["gene_symbol", "title", "summary"]
    all_cols  = base_cols + (extra_cols or []) + ["payload_str"]

    rows = []
    for r in records:
        row = (
            r.get("gene_symbol", ""),
            r.get("title", ""),
            r.get("summary", ""),
        )
        for c in (extra_cols or []):
            row += (r.get(c, ""),)
        row += (json.dumps(r.get("payload", {}), default=str),)
        rows.append(row)

    schema_fields = [
        StructField(c, StringType(), True)
        for c in all_cols
    ]
    df_raw = spark.createDataFrame(rows, StructType(schema_fields))

    select_exprs = [F.col(c) for c in all_cols[:-1]]  # all non-payload cols
    select_exprs.append(F.expr("parse_json(payload_str)").alias("payload"))

    df_final = df_raw.select(*select_exprs)
    print(f"  Rows to write : {df_final.count():,}")
    (
        df_final.write
        .format("delta")
        .mode("overwrite")
        .option("overwriteSchema", "false")
        .saveAsTable(table)
    )
    print(f"  ✅ Written to {table}")


print("Writing neuroplex_uniprot ...")
_write_records(uniprot_records, TABLE_UNIPROT)

# COMMAND ----------

# DBTITLE 1,Write neuroplex_kegg to Delta
print("Writing neuroplex_kegg ...")
_write_records(kegg_records, TABLE_KEGG)

# COMMAND ----------

# DBTITLE 1,Validation — row counts + PSEN1 spot-check
# MAGIC %sql
# MAGIC -- UniProt: should have 1 row per gene (canonical SwissProt entry)
# MAGIC SELECT gene_symbol,
# MAGIC        payload:uniprot_acc::STRING                  AS accession,
# MAGIC        payload:protein_name::STRING                 AS protein,
# MAGIC        LEFT(payload:function::STRING, 120)          AS function_snippet,
# MAGIC        payload:subcellular_location::STRING         AS location
# MAGIC FROM dhbl_discovery_us_val.genesis_schema.neuroplex_uniprot
# MAGIC ORDER BY gene_symbol

# COMMAND ----------

# DBTITLE 1,KEGG validation — PSEN1 pathway membership
# MAGIC %sql
# MAGIC -- KEGG: PSEN1 should appear in Alzheimer disease (hsa05010) and Notch signaling
# MAGIC SELECT gene_symbol,
# MAGIC        payload:pathway_id::STRING    AS pathway_id,
# MAGIC        payload:pathway_name::STRING  AS pathway_name
# MAGIC FROM dhbl_discovery_us_val.genesis_schema.neuroplex_kegg
# MAGIC WHERE UPPER(gene_symbol) = 'PSEN1'
# MAGIC ORDER BY pathway_id

# COMMAND ----------

# DBTITLE 1,GBA1 patch — append UniProt + KEGG rows
# Run cells 1-3 first if kernel was restarted.
# Appends GBA1 rows to neuroplex_uniprot and neuroplex_kegg.

import json as _json
from pyspark.sql import functions as _F
from pyspark.sql.types import StructType as _ST, StructField as _SF, StringType as _Str

def _append_gba1(records, table, str_cols):
    if not records:
        print(f"  ⚠️  No GBA1 records → {table}")
        return
    rows = [tuple(r.get(c, "") for c in str_cols)
            + (_json.dumps(r.get("payload", {}), default=str),) for r in records]
    schema = _ST([_SF(c, _Str(), True) for c in str_cols] + [_SF("payload_str", _Str(), True)])
    df = spark.createDataFrame(rows, schema).select(
        *[_F.col(c) for c in str_cols],
        _F.expr("parse_json(payload_str)").alias("payload")
    )
    n = df.count()
    df.write.format("delta").mode("append").saveAsTable(table)
    print(f"  ✅ Appended {n} GBA1 rows → {table}")


print("GBA1 → UniProt ...")
u = fetch_uniprot("GBA1")
if u: print(f"  {u['payload']['uniprot_acc']}  {u['payload']['protein_name']}")
_append_gba1([u] if u else [], TABLE_UNIPROT, ["gene_symbol", "title", "summary"])

time.sleep(REQUEST_DELAY)

print("\nGBA1 → KEGG ...")
kr = fetch_kegg_pathways("GBA1")
print(f"  {len(kr)} pathways")
_append_gba1(kr, TABLE_KEGG, ["gene_symbol", "title", "summary"])

# COMMAND ----------


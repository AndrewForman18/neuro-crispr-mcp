# Databricks notebook source
# /// script
# [tool.databricks.environment]
# environment_version = "5"
# ///
# DBTITLE 1,Config — environment, gene panel, API endpoint
import os, sys, json, time, logging
import requests
from pyspark.sql import functions as F
from pyspark.sql.types import StructType, StructField, StringType

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("opentargets_etl")

os.environ["NEUROPLEX_ENV"] = "val"
sys.path.insert(0, "/Workspace/Users/andrew_forman@eisai.com/neuro-crispr-mcp")
for mod in [k for k in sys.modules if k.startswith("config")]:
    del sys.modules[mod]
from config.neuroplex_config import load_config

CFG = load_config()
CAT = CFG.catalog
SCH = CFG.schema
TABLE = f"{CAT}.{SCH}.neuroplex_opentargets"

print(f"✅ Environment : {CFG.environment}")
print(f"   Catalog     : {CAT}.{SCH}")
print(f"   Target table: {TABLE}")

GENE_PANEL = [
    "PSEN1", "PSEN2", "APP",   "APOE",  "TREM2",  "BIN1",
    "CLU",   "CR1",   "PICALM","ABCA7", "SORL1",
    "SNCA",  "LRRK2", "PINK1", "PRKN",  "PARK7",  "GBA1", "VPS35",
    "SOD1",  "TARDBP","FUS",   "C9orf72","TBK1",  "GRN",
    "MAPT",  "HTT",
    "HCRT",  "HCRTR1","HCRTR2",
    "MECP2", "FMR1",
]

OT_API       = "https://api.platform.opentargets.org/api/v4/graphql"
REQUEST_DELAY = 1.0

print(f"   Gene panel  : {len(GENE_PANEL)} genes")

# COMMAND ----------

# DBTITLE 1,OpenTargets GraphQL helpers
# OpenTargets Platform v4 GraphQL API (public, no auth required)
# Pipeline per gene:
#   1. Search gene symbol → resolve Ensembl target ID (ENSG...)
#   2. Query associatedDiseases for that target → disease name + OT score
# payload shape matches neuroplex_query.py selectors:
#   payload:score::DOUBLE
#   payload:disease.therapeuticAreas[0].name::STRING

# ── GraphQL queries ──────────────────────────────────────────────────────────────
SEARCH_TARGET_QUERY = """
query SearchTarget($q: String!) {
  search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 1}) {
    hits {
      id
      name
    }
  }
}
"""

ASSOCIATIONS_QUERY = """
query GeneAssociations($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    approvedSymbol
    associatedDiseases(
      enableIndirect: true
      page: {index: 0, size: 100}
    ) {
      rows {
        disease {
          id
          name
          therapeuticAreas {
            id
            name
          }
        }
        score
        datatypeScores {
          id
          score
        }
      }
    }
  }
}
"""


def _ot_gql(query: str, variables: dict, retries: int = 3) -> dict | None:
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    for attempt in range(retries):
        try:
            r = requests.post(OT_API, json={"query": query, "variables": variables},
                              headers=headers, timeout=30)
            if r.status_code == 200:
                body = r.json()
                if "errors" in body:
                    log.warning("OT GraphQL errors: %s", body["errors"])
                    return None
                return body.get("data")
            elif r.status_code == 429:
                wait = 2 ** (attempt + 1)
                log.warning("OT rate limited, waiting %ss", wait)
                time.sleep(wait)
            else:
                log.warning("OT HTTP %s", r.status_code)
                return None
        except Exception as exc:
            log.warning("OT request error (attempt %d): %s", attempt + 1, exc)
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
    return None


def resolve_ensembl_id(gene_symbol: str) -> str | None:
    """Resolve gene symbol to an Ensembl target ID via OpenTargets search."""
    data = _ot_gql(SEARCH_TARGET_QUERY, {"q": gene_symbol})
    if not data:
        return None
    hits = (data.get("search") or {}).get("hits") or []
    if not hits:
        return None
    # Verify the top hit name matches (OT search is substring-based)
    top = hits[0]
    if top.get("name", "").upper() == gene_symbol.upper():
        return top["id"]   # ENSG...
    return None


def fetch_ot_associations(gene_symbol: str) -> list[dict]:
    """
    Fetch gene-disease associations from OpenTargets for a single gene.
    Returns a list of row dicts ready for neuroplex_opentargets.
    """
    ensembl_id = resolve_ensembl_id(gene_symbol)
    if not ensembl_id:
        log.warning("Could not resolve Ensembl ID for %s", gene_symbol)
        return []

    time.sleep(REQUEST_DELAY)

    data = _ot_gql(ASSOCIATIONS_QUERY, {"ensemblId": ensembl_id})
    if not data or not data.get("target"):
        return []

    rows = ((data["target"].get("associatedDiseases") or {}).get("rows") or [])
    records = []
    for row in rows:
        disease_obj  = row.get("disease") or {}
        disease_name = disease_obj.get("name", "")
        score        = row.get("score")
        datatype_scores = {s["id"]: s["score"] for s in (row.get("datatypeScores") or [])}

        records.append({
            "gene_symbol": gene_symbol,
            "disease": disease_name,
            "payload": {
                "ensembl_id": ensembl_id,
                "score": score,
                "disease": {
                    "id":   disease_obj.get("id"),
                    "name": disease_name,
                    "therapeuticAreas": disease_obj.get("therapeuticAreas") or [],
                },
                "datatype_scores": datatype_scores,
            },
        })
    return records


print("✅ OpenTargets helpers defined")

# COMMAND ----------

# DBTITLE 1,API connectivity test — PSEN1 probe
print("Resolving PSEN1 Ensembl ID ...")
eid = resolve_ensembl_id("PSEN1")
print(f"  Ensembl ID : {eid}")

if eid:
    time.sleep(REQUEST_DELAY)
    print("\nFetching PSEN1 disease associations ...")
    recs = fetch_ot_associations("PSEN1")
    print(f"  {len(recs)} associations found")
    if recs:
        for r in sorted(recs, key=lambda x: x['payload']['score'], reverse=True)[:5]:
            area = (r['payload']['disease']['therapeuticAreas'] or [{}])[0].get('name', '')
            print(f"  {r['payload']['score']:.3f}  {r['disease'][:45]:45s}  [{area}]")
        print("  ✅ API OK — ready to run full panel")
else:
    print("  ⚠️  Could not resolve PSEN1 — check OT search endpoint")

# COMMAND ----------

# DBTITLE 1,Fetch OT associations for all panel genes
all_records = []
errors      = []

for i, gene in enumerate(GENE_PANEL):
    print(f"[{i+1:02d}/{len(GENE_PANEL)}] {gene} ...", end=" ")
    try:
        recs = fetch_ot_associations(gene)
        all_records.extend(recs)
        if recs:
            top = max(recs, key=lambda r: r['payload']['score'] or 0)
            print(f"✅  {len(recs)} diseases  (top: {top['disease'][:40]!r} score={top['payload']['score']:.3f})")
        else:
            print("⚠️  0 associations")
    except Exception as exc:
        errors.append((gene, str(exc)))
        print(f"❌  {exc}")
    time.sleep(REQUEST_DELAY * 2)  # OT API: 2s between genes (search + assoc queries)

print(f"\n———")
print(f"Total records : {len(all_records):,}")
print(f"Genes covered : {len(set(r['gene_symbol'] for r in all_records))}/{len(GENE_PANEL)}")
if errors:
    print(f"Errors ({len(errors)}): {errors}")

# COMMAND ----------

# DBTITLE 1,Write neuroplex_opentargets to Delta
from pyspark.sql import functions as F
from pyspark.sql.types import StructType, StructField, StringType

if not all_records:
    raise RuntimeError("No records fetched — check API connectivity before writing.")

# Schema: gene_symbol | disease | payload_str  →  cast to VARIANT
rows = [
    (
        r["gene_symbol"],
        r["disease"],
        json.dumps(r["payload"], default=str),
    )
    for r in all_records
]

df_raw = spark.createDataFrame(
    rows,
    StructType([
        StructField("gene_symbol",  StringType(), True),
        StructField("disease",       StringType(), True),
        StructField("payload_str",   StringType(), True),
    ])
)

df_final = df_raw.select(
    F.col("gene_symbol"),
    F.col("disease"),
    F.expr("parse_json(payload_str)").alias("payload"),
)

print(f"Rows to write : {df_final.count():,}")

(
    df_final.write
    .format("delta")
    .mode("overwrite")
    .option("overwriteSchema", "false")
    .saveAsTable(TABLE)
)

print(f"✅ Written to {TABLE}")

# COMMAND ----------

# DBTITLE 1,Validation — PSEN1 disease associations
# MAGIC %sql
# MAGIC -- PSEN1 should have strong Alzheimer associations (score close to 1.0)
# MAGIC SELECT
# MAGIC     gene_symbol,
# MAGIC     disease,
# MAGIC     payload:score::DOUBLE                                   AS ot_score,
# MAGIC     payload:disease.therapeuticAreas[0].name::STRING        AS therapeutic_area
# MAGIC FROM dhbl_discovery_us_val.genesis_schema.neuroplex_opentargets
# MAGIC WHERE UPPER(gene_symbol) = 'PSEN1'
# MAGIC ORDER BY ot_score DESC
# MAGIC LIMIT 20

# COMMAND ----------

# DBTITLE 1,Validation — panel-wide disease landscape
# MAGIC %sql
# MAGIC -- Top-scoring gene-disease pairs across the full NeuroPlex panel
# MAGIC SELECT
# MAGIC     gene_symbol,
# MAGIC     disease,
# MAGIC     payload:score::DOUBLE                             AS ot_score,
# MAGIC     payload:disease.therapeuticAreas[0].name::STRING  AS therapeutic_area
# MAGIC FROM dhbl_discovery_us_val.genesis_schema.neuroplex_opentargets
# MAGIC WHERE payload:score::DOUBLE >= 0.5
# MAGIC ORDER BY ot_score DESC
# MAGIC LIMIT 30

# COMMAND ----------

# DBTITLE 1,GBA1 patch — append OT associations
# Run cells 1-2 first if kernel was restarted.
# Appends GBA1 disease association rows to neuroplex_opentargets.

import json as _json
from pyspark.sql import functions as _F
from pyspark.sql.types import StructType as _ST, StructField as _SF, StringType as _Str

print("GBA1 → OpenTargets ...")
gba1_recs = fetch_ot_associations("GBA1")
if gba1_recs:
    top = max(gba1_recs, key=lambda r: r['payload']['score'] or 0)
    print(f"  {len(gba1_recs)} associations  (top: {top['disease']!r} score={top['payload']['score']:.3f})")
else:
    print("  ⚠️  0 associations")

if gba1_recs:
    rows = [(r["gene_symbol"], r["disease"],
             _json.dumps(r["payload"], default=str)) for r in gba1_recs]
    schema = _ST([_SF("gene_symbol", _Str(), True), _SF("disease", _Str(), True),
                  _SF("payload_str", _Str(), True)])
    df = spark.createDataFrame(rows, schema).select(
        _F.col("gene_symbol"), _F.col("disease"),
        _F.expr("parse_json(payload_str)").alias("payload")
    )
    n = df.count()
    df.write.format("delta").mode("append").saveAsTable(TABLE)
    print(f"  ✅ Appended {n} GBA1 rows → {TABLE}")

# COMMAND ----------


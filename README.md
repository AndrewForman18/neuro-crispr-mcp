# NeuroPlex

NeuroPlex is a Databricks App for exploring cross-species CRISPR perturbation, expression, pharmacology, genetics, protein, and pathway datasets through a Streamlit UI backed by Unity Catalog tables and Databricks model serving.

## Repository contents

* `app.py` — Streamlit app UI and orchestration
* `app.yaml` — Databricks App manifest for the active environment
* `app.yaml.template` — manifest template reference
* `config/neuroplex_env.yml` — environment-specific settings
* `config/neuroplex_config.py` — runtime config loader and `app.yaml` renderer
* `neuro_mcp_server/` — SQL-backed query helpers and chart generation
* `ingestion/` — source registry, ingestors, and ingestion notebooks/tasks
* `scripts/render_app_yaml.py` — renders `app.yaml` for a chosen environment

## Multi-environment packaging

Environment-specific values were abstracted into `config/neuroplex_env.yml`.

Current runtime config reads:

* `catalog`
* `schema`
* `sql_warehouse_id`
* `serving_endpoint`
* `table_prefix`
* `query_log_table`
* app resource bindings for SQL Warehouse and serving endpoint

You can also override key values at deploy/runtime with environment variables:

* `NEUROPLEX_ENV`
* `NEUROPLEX_CATALOG`
* `NEUROPLEX_SCHEMA`
* `NEUROPLEX_SQL_WAREHOUSE_ID`
* `NEUROPLEX_TABLE_PREFIX`
* `NEUROPLEX_QUERY_LOG_TABLE`
* `DATABRICKS_SERVING_ENDPOINT`
* `DATABRICKS_HOST`

## Set up a new environment

1. Edit `config/neuroplex_env.yml` and populate the target environment section.
2. Render the Databricks App manifest:

```bash
PYTHONPATH=. python scripts/render_app_yaml.py --env staging --output app.yaml
```

3. Verify the rendered `app.yaml` has the correct SQL warehouse and serving endpoint bindings.
4. Deploy from a Databricks workspace folder containing this repository.

## Databricks Apps deploy sequence

After updating source files, use this sequence:

```bash
databricks apps get neuro-crispr-mcp --output JSON
```

Then branch on state:

* `RUNNING` with no `pending_deployment`: deploy
* `STOPPED`: start first, then deploy
* `STARTING` or `STOPPING`: poll `apps get` until stable

Deploy example:

```bash
databricks apps deploy neuro-crispr-mcp --source-code-path /Workspace/Users/<user>/neuro-crispr-mcp --output JSON
```

## Notes

* Exported notebooks are stored as `.ipynb` files for Git compatibility.
* Some notebook output cells still contain historical dev-environment text from prior runs; the runtime and primary configuration paths are now environment-driven.
* `requirements.txt` now includes `PyYAML` for config loading.

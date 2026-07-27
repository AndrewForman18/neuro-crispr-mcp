from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

CONFIG_PATH = Path(__file__).with_name("neuroplex_env.yml")

@dataclass(frozen=True)
class NeuroPlexConfig:
    environment: str
    workspace_host: str
    catalog: str
    schema: str
    sql_warehouse_id: str
    serving_endpoint: str
    table_prefix: str = "neuroplex"
    query_log_table: str = "neuroplex_query_log"
    app_resources: dict[str, Any] | None = None

    def fqn(self, table_name: str) -> str:
        return f"{self.catalog}.{self.schema}.{table_name}"

    def prefixed_table(self, suffix: str) -> str:
        return f"{self.table_prefix}_{suffix}"

    def prefixed_fqn(self, suffix: str) -> str:
        return self.fqn(self.prefixed_table(suffix))

    @property
    def query_log_fqn(self) -> str:
        return self.fqn(self.query_log_table)

    def render_app_yaml(self) -> str:
        wh = (self.app_resources or {}).get("sql_warehouse", {})
        se = (self.app_resources or {}).get("serving_endpoint", {})
        return (
            "command:\n"
            "  - streamlit\n"
            "  - run\n"
            "  - app.py\n"
            "  - --server.port\n"
            "  - \"8000\"\n"
            "env:\n"
            "  - name: STREAMLIT_GATHER_USAGE_STATS\n"
            "    value: \"false\"\n"
            f"  - name: NEUROPLEX_ENV\n    value: {self.environment}\n"
            "  - name: DATABRICKS_SQL_WAREHOUSE_HTTP_PATH\n"
            f"    value: /sql/1.0/warehouses/{wh.get('id', self.sql_warehouse_id)}\n"
            f"  - name: DATABRICKS_HOST\n    value: {self.workspace_host}\n"
            "  - name: DATABRICKS_SERVING_ENDPOINT\n"
            f"    value: {se.get('name', self.serving_endpoint)}\n"
            "resources:\n"
            "  - name: sql-warehouse\n"
            "    sql_warehouse:\n"
            f"      id: {wh.get('id', self.sql_warehouse_id)}\n"
            f"      permission: {wh.get('permission', 'CAN_USE')}\n"
            "  - name: serving-endpoint\n"
            "    serving_endpoint:\n"
            f"      name: {se.get('name', self.serving_endpoint)}\n"
            f"      permission: {se.get('permission', 'CAN_QUERY')}\n"
        )


def _load_raw() -> dict[str, Any]:
    with CONFIG_PATH.open() as f:
        return yaml.safe_load(f)


def load_config(env: str | None = None) -> NeuroPlexConfig:
    raw = _load_raw()
    active_env = env or os.environ.get("NEUROPLEX_ENV") or raw.get("active_environment", "dev")
    data = dict(raw["environments"][active_env])
    data["workspace_host"] = os.environ.get("DATABRICKS_HOST", data["workspace_host"])
    data["catalog"] = os.environ.get("NEUROPLEX_CATALOG", data["catalog"])
    data["schema"] = os.environ.get("NEUROPLEX_SCHEMA", data["schema"])
    data["sql_warehouse_id"] = os.environ.get("NEUROPLEX_SQL_WAREHOUSE_ID", data["sql_warehouse_id"])
    data["serving_endpoint"] = os.environ.get("DATABRICKS_SERVING_ENDPOINT", data["serving_endpoint"])
    data["table_prefix"] = os.environ.get("NEUROPLEX_TABLE_PREFIX", data.get("table_prefix", "neuroplex"))
    data["query_log_table"] = os.environ.get("NEUROPLEX_QUERY_LOG_TABLE", data.get("query_log_table", "neuroplex_query_log"))
    return NeuroPlexConfig(environment=active_env, **data)

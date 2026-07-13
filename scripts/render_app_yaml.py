from __future__ import annotations

import argparse
from pathlib import Path

from config.neuroplex_config import load_config

parser = argparse.ArgumentParser(description="Render app.yaml for a target NeuroPlex environment")
parser.add_argument("--env", default=None)
parser.add_argument("--output", default="app.yaml")
args = parser.parse_args()

cfg = load_config(args.env)
Path(args.output).write_text(cfg.render_app_yaml())
print(f"Rendered {args.output} for environment={cfg.environment}")

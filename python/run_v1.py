from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from hydroiso.config import ConfigurationError, apply_overrides, load_config
from hydroiso.runner import execute, execute_grid_convergence


def main() -> int:
    parser = argparse.ArgumentParser(description="Run hydroisomerization pellet model V1")
    parser.add_argument("config", type=Path, help="User YAML configuration")
    parser.add_argument("--set", dest="overrides", action="append", default=[], metavar="PATH=VALUE")
    parser.add_argument("--build-dir", type=Path, default=Path("build"))
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--grid-convergence", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    try:
        data = load_config(args.config.resolve())
        data = apply_overrides(data, args.overrides)
        run_directory = execute(root, data, (root / args.build_dir).resolve(), args.skip_build)
        if args.grid_convergence or data["validation"]["grid_convergence"]["enabled"]:
            table = execute_grid_convergence(root, data, (root / args.build_dir).resolve(), run_directory)
            print(f"Grid convergence: {table}")
    except (ConfigurationError, RuntimeError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(f"Completed: {run_directory}")
    if data["run_mode"] == "demo":
        print("WARNING: RESULTS ARE NOT EXPERIMENTALLY VALIDATED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

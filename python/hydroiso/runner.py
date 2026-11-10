from __future__ import annotations

import json
import subprocess
import time
import copy
import csv
import hashlib
from pathlib import Path
from typing import Any

from .config import dump_normalized, normalized_hash, write_fortran_namelist
from .plots import create_plots


def build_project(root: Path, build_directory: Path) -> Path:
    cache = build_directory / "CMakeCache.txt"
    if not cache.exists():
        _run([
            "cmake", "-S", str(root), "-B", str(build_directory), "-G", "Ninja",
            "-DCMAKE_BUILD_TYPE=Debug", "-DCMAKE_PREFIX_PATH=C:/msys64/ucrt64",
        ], root)
    _run(["cmake", "--build", str(build_directory)], root)
    executable = build_directory / "hydroiso_v1.exe"
    if not executable.exists():
        raise RuntimeError(f"Fortran executable was not produced: {executable}")
    return executable


def execute(
    root: Path, data: dict[str, Any], build_directory: Path, skip_build: bool = False,
    make_plots: bool = True,
) -> Path:
    digest = normalized_hash(data)
    source_digest = _source_hash(root)
    run_directory = root / data["output"]["root_directory"] / f"{data['run_name']}_{digest}_{source_digest}"
    if run_directory.exists() and not data["output"].get("overwrite_existing_run", False):
        if (run_directory / "run_complete.json").exists():
            return run_directory
        raise RuntimeError(f"Incomplete run directory already exists: {run_directory}")
    run_directory.mkdir(parents=True, exist_ok=True)
    dump_normalized(data, run_directory / "normalized_config.yaml")
    write_fortran_namelist(data, run_directory / "fortran_config.nml")
    executable = build_directory / "hydroiso_v1.exe" if skip_build else build_project(root, build_directory)
    started = time.time()
    completed = subprocess.run(
        [str(executable), str(run_directory / "fortran_config.nml"), str(run_directory)],
        cwd=root, text=True, encoding="utf-8", errors="replace", capture_output=True,
    )
    (run_directory / "stdout.log").write_text(completed.stdout, encoding="utf-8")
    (run_directory / "stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"Fortran calculation failed with code {completed.returncode}. "
            f"See {run_directory / 'stderr.log'}"
        )
    plot_paths = create_plots(run_directory, data) if make_plots else []
    metadata = {
        "schema_version": data["schema_version"],
        "model_version": data["model_version"],
        "run_mode": data["run_mode"],
        "config_hash": digest,
        "source_hash": source_digest,
        "elapsed_seconds": time.time()-started,
        "warning": "RESULTS ARE NOT EXPERIMENTALLY VALIDATED" if data["run_mode"] == "demo" else None,
        "plots": [path.name for path in plot_paths],
    }
    (run_directory / "run_complete.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return run_directory


def _source_hash(root: Path) -> str:
    digest = hashlib.sha256()
    paths = [root / "CMakeLists.txt"]
    paths.extend(sorted((root / "src").rglob("*.f90")))
    paths.extend(sorted((root / "python" / "hydroiso").glob("*.py")))
    for path in paths:
        digest.update(path.relative_to(root).as_posix().encode("utf-8"))
        digest.update(path.read_bytes())
    return digest.hexdigest()[:10]


def execute_grid_convergence(
    root: Path, data: dict[str, Any], build_directory: Path, parent_run: Path,
) -> Path:
    grids = [int(value) for value in data["validation"]["grid_convergence"]["radial_cells"]]
    grids = sorted(set(grids))
    rows: list[dict[str, float | int | str]] = []
    previous_conversion: float | None = None
    tolerance = float(data["validation"]["grid_convergence"]["relative_tolerance"])
    for grid in grids:
        case = copy.deepcopy(data)
        case["run_name"] = f"{data['run_name']}_grid_{grid}"
        case["numerics"]["spatial"]["radial_cells"] = grid
        case["numerics"]["time"]["output_interval"]["value"] = case["numerics"]["time"]["end_time"]["value"]
        case["validation"]["grid_convergence"]["enabled"] = False
        run_dir = execute(root, case, build_directory, skip_build=True, make_plots=False)
        conversion = _final_conversion(run_dir)
        error = 0.0 if previous_conversion is None else abs(conversion-previous_conversion)/max(abs(conversion), 1.0e-15)
        rows.append({
            "radial_cells": grid,
            "conversion": conversion,
            "relative_change": error,
            "converged": previous_conversion is not None and error <= tolerance,
            "run_directory": run_dir.name,
        })
        previous_conversion = conversion
    path = parent_run / "grid_convergence.csv"
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    _plot_grid_convergence(path, parent_run / "plots" / "grid_convergence.png")
    return path


def _final_conversion(run_directory: Path) -> float:
    with (run_directory / "averages.csv").open("r", encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    initial = float(rows[0]["n-C16H34"])
    final = float(rows[-1]["n-C16H34"])
    return 1.0-final/initial


def _plot_grid_convergence(csv_path: Path, output_path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with csv_path.open("r", encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    grids = [int(row["radial_cells"]) for row in rows]
    changes = [float(row["relative_change"]) for row in rows]
    fig, axis = plt.subplots(figsize=(7, 4.5))
    axis.semilogy(grids[1:], changes[1:], marker="o")
    axis.set(xlabel="radial cells", ylabel="relative change in final conversion")
    axis.grid(True, which="both", alpha=0.3)
    fig.tight_layout(); fig.savefig(output_path, dpi=160); plt.close(fig)


def _run(command: list[str], cwd: Path) -> None:
    completed = subprocess.run(command, cwd=cwd, text=True, encoding="utf-8", errors="replace")
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed ({completed.returncode}): {' '.join(command)}")

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.linalg import solve_banded
from scipy.special import ive


def create_plots(run_directory: Path, config: dict) -> list[Path]:
    plot_dir = run_directory / "plots"
    plot_dir.mkdir(exist_ok=True)
    profiles = pd.read_csv(run_directory / "profiles.csv")
    averages = pd.read_csv(run_directory / "averages.csv")
    rates = pd.read_csv(run_directory / "reaction_rates.csv")
    selected = config["plots"].get("selected_species", ["n-C16H34", "H2"])
    species_columns = list(profiles.columns[3:])
    selected_columns = _resolve_species(selected, species_columns)
    outputs: list[Path] = []

    selected_times = config["plots"]["selected_times"]["value"]
    fig, axes = plt.subplots(len(selected_columns), 1, figsize=(7, 3.2*len(selected_columns)), squeeze=False)
    available_times = profiles["time_s"].unique()
    for axis, column in zip(axes[:, 0], selected_columns):
        for requested in selected_times:
            time = available_times[np.argmin(abs(available_times-float(requested)))]
            frame = profiles[np.isclose(profiles["time_s"], time)]
            axis.plot(frame["r_m"]*1e3, frame[column], label=f"t={time:g} s")
        axis.set_ylabel(f"{column}, mol/m³")
        axis.grid(True, alpha=0.3)
        axis.legend()
    axes[-1, 0].set_xlabel("r, mm")
    outputs.append(_save(fig, plot_dir / "concentration_profiles.png"))

    center = profiles.loc[profiles.groupby("time_s")["cell"].idxmin()]
    surface = profiles.loc[profiles.groupby("time_s")["cell"].idxmax()]
    fig, axes = plt.subplots(len(selected_columns), 1, figsize=(7, 3.2*len(selected_columns)), squeeze=False)
    for axis, column in zip(axes[:, 0], selected_columns):
        axis.plot(center["time_s"], center[column], label="center")
        axis.plot(surface["time_s"], surface[column], label="surface cell")
        axis.set_ylabel(f"{column}, mol/m³")
        axis.grid(True, alpha=0.3)
        axis.legend()
    axes[-1, 0].set_xlabel("time, s")
    for axis in axes[:, 0]:
        axis.set_xscale("symlog", linthresh=1.0)
    outputs.append(_save(fig, plot_dir / "center_surface_vs_time.png"))

    fig, axes = plt.subplots(len(selected_columns), 1, figsize=(7, 3.2*len(selected_columns)), squeeze=False)
    for axis, column in zip(axes[:, 0], selected_columns):
        axis.plot(averages["time_s"], averages[column], label=column)
        axis.set_ylabel("mol/m³")
        axis.set_xscale("symlog", linthresh=1.0)
        axis.grid(True, alpha=0.3); axis.legend()
    axes[-1, 0].set_xlabel("time, s")
    outputs.append(_save(fig, plot_dir / "average_concentrations.png"))

    fig, axis = plt.subplots(figsize=(8, 5))
    for column in [f"R{i}" for i in range(1, 16)]:
        axis.plot(rates["time_s"], rates[column], label=column)
    axis.set(xlabel="time, s", ylabel="volume-average reaction rate, mol/(m³ s)")
    axis.set_xscale("symlog", linthresh=1.0)
    axis.set_yscale("symlog", linthresh=1.0e-5)
    axis.grid(True, alpha=0.3); axis.legend(ncol=3, fontsize=8)
    outputs.append(_save(fig, plot_dir / "reaction_rates.png"))

    fig, axis = plt.subplots(figsize=(8, 5))
    for column in [f"eta_R{i}" for i in range(1, 16)]:
        values = pd.to_numeric(rates[column], errors="coerce").replace([np.inf, -np.inf], np.nan)
        if np.any(np.abs(values.fillna(0.0).to_numpy()) > 0.0):
            axis.plot(rates["time_s"], values, label=column)
    axis.set(xlabel="time, s", ylabel="numerical effectiveness factor")
    axis.set_xscale("symlog", linthresh=1.0)
    axis.grid(True, alpha=0.3); axis.legend(ncol=3, fontsize=8)
    outputs.append(_save(fig, plot_dir / "effectiveness_factors.png"))

    initial_c1 = float(averages.iloc[0]["n-C16H34"])
    conversion = 1.0 - averages["n-C16H34"]/initial_c1
    pd.DataFrame({"time_s": averages["time_s"], "conversion_C16": conversion}).to_csv(
        run_directory / "conversion.csv", index=False
    )
    fig, axis = plt.subplots(figsize=(7, 4.5))
    axis.plot(averages["time_s"], conversion)
    axis.set(xlabel="time, s", ylabel="pellet-average n-C16 conversion")
    axis.set_xscale("symlog", linthresh=1.0)
    axis.grid(True, alpha=0.3)
    outputs.append(_save(fig, plot_dir / "conversion.png"))
    if config["plots"].get("thiele_test_diagnostic", False):
        outputs.append(_thiele_diagnostic(run_directory, plot_dir))
    return outputs


def _thiele_diagnostic(run_directory: Path, plot_dir: Path) -> Path:
    phi = np.logspace(-3, 2, 80)
    analytical = 2.0*ive(1, phi)/(phi*ive(0, phi))
    numerical = np.array([_first_order_cylinder_eta(value, 200) for value in phi])
    relative_error = np.abs(numerical-analytical)/np.maximum(np.abs(analytical), 1.0e-15)
    pd.DataFrame({
        "phi": phi,
        "eta_analytical": analytical,
        "eta_fvm_200": numerical,
        "relative_error": relative_error,
    }).to_csv(run_directory / "thiele_diagnostic.csv", index=False)
    fig, axes = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
    axes[0].loglog(phi, analytical, label="analytical cylinder")
    axes[0].loglog(phi, numerical, "--", label="FVM, Nr=200")
    axes[0].set_ylabel("effectiveness factor")
    axes[0].grid(True, which="both", alpha=0.3); axes[0].legend()
    axes[1].loglog(phi, relative_error)
    axes[1].set(xlabel="Thiele modulus", ylabel="relative error")
    axes[1].grid(True, which="both", alpha=0.3)
    return _save(fig, plot_dir / "thiele_diagnostic.png")


def _first_order_cylinder_eta(phi: float, cells: int) -> float:
    dr = 1.0/cells
    faces = np.arange(cells+1, dtype=float)*dr
    areas = 2.0*np.pi*faces
    volumes = np.pi*(faces[1:]**2-faces[:-1]**2)
    banded = np.zeros((3, cells))
    rhs = np.zeros(cells)
    for j in range(cells):
        g_left = 0.0 if j == 0 else areas[j]/dr
        g_right = areas[j+1]/dr if j < cells-1 else areas[j+1]/(0.5*dr)
        banded[1, j] = -(g_left+g_right+phi*phi*volumes[j])
        if j > 0:
            banded[2, j-1] = g_left
        if j < cells-1:
            banded[0, j+1] = g_right
        else:
            rhs[j] = -g_right
    concentration = solve_banded((1, 1), banded, rhs)
    return float(np.sum(concentration*volumes)/np.sum(volumes))


def _resolve_species(requested: list[str], columns: list[str]) -> list[str]:
    result = []
    for item in requested:
        if item in columns:
            result.append(item)
        elif item.startswith("C") and item[1:].isdigit():
            index = int(item[1:]) - 1
            if 0 <= index < len(columns):
                result.append(columns[index])
    return list(dict.fromkeys(result))


def _save(fig: plt.Figure, path: Path) -> Path:
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path

from __future__ import annotations

import copy
import hashlib
import json
import re
import math
from pathlib import Path
from typing import Any

import yaml


class ConfigurationError(ValueError):
    """Raised when a user configuration is incomplete or inconsistent."""


def load_config(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        data = yaml.safe_load(stream)
    if not isinstance(data, dict):
        raise ConfigurationError("The YAML root must be a mapping.")
    data = normalize_derived(data)
    validate_config(data)
    return data


def apply_overrides(data: dict[str, Any], overrides: list[str]) -> dict[str, Any]:
    result = copy.deepcopy(data)
    for expression in overrides:
        if "=" not in expression:
            raise ConfigurationError(f"Override must be path=value, got: {expression}")
        path, raw_value = expression.split("=", 1)
        keys = path.split(".")
        node: Any = result
        for key in keys[:-1]:
            if not isinstance(node, dict) or key not in node:
                raise ConfigurationError(f"Unknown override path: {path}")
            node = node[key]
        if not isinstance(node, dict) or keys[-1] not in node:
            raise ConfigurationError(f"Unknown override path: {path}")
        node[keys[-1]] = yaml.safe_load(raw_value)
    result = normalize_derived(result)
    validate_config(result)
    return result


def normalize_derived(data: dict[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(data)
    pellet = result["pellet"]
    vpore = float(pellet["micropore_volume"]["value"]) + float(pellet["mesopore_volume"]["value"])
    pellet["total_pore_volume"]["value"] = vpore
    if pellet["porosity"]["mode"] == "from_skeletal_density_and_pore_volume":
        rho_s = float(pellet["skeletal_density"]["value"])
        pellet["porosity"]["value"] = vpore/(vpore+1.0/rho_s)
        pellet["particle_density"]["value"] = 1.0/(vpore+1.0/rho_s)
    if pellet["pore_radius"]["mode"] == "calculated_from_total_pore_volume_and_bet":
        pellet["pore_radius"]["value"] = 2.0*vpore/float(pellet["bet_surface_area"]["value"])
    kinetics = result["kinetics"]
    normalization = kinetics["demo_embedding"]["normalization_concentration"]
    if normalization["mode"] == "ideal_gas_total_at_bulk_conditions":
        pressure = float(result["conditions"]["pressure"]["value"])
        temperature = float(result["conditions"]["temperature"]["value"])
        normalization["value"] = pressure/(8.314462618*temperature)
    ctotal = float(normalization["value"])
    result["bulk"]["total_concentration"]["value"] = ctotal
    if result["bulk"].get("composition_input") == "mole_fraction":
        fractions = result["bulk"]["mole_fractions"]
        for i in range(1, 16):
            result["bulk"]["concentrations"][f"C{i}"]["value"] = ctotal*float(fractions[f"C{i}"])
    kinetics["acid_sites"]["estimated_bronsted_concentration"]["value"] = (
        float(kinetics["acid_sites"]["bronsted_inventory"]["value"])
        * float(pellet["particle_density"]["value"])
    )
    return result


def validate_config(data: dict[str, Any]) -> None:
    required = [
        "schema_version", "model_version", "run_mode", "conditions", "pellet",
        "species", "kinetics", "diffusion", "bulk", "initial_conditions",
        "numerics", "output", "plots",
    ]
    missing = [key for key in required if key not in data]
    if missing:
        raise ConfigurationError(f"Missing top-level sections: {', '.join(missing)}")
    if data["model_version"] != "V1":
        raise ConfigurationError("This runner implements model_version: V1 only.")
    if data["run_mode"] not in {"demo", "scientific"}:
        raise ConfigurationError("run_mode must be demo or scientific.")
    if data["run_mode"] == "scientific" and _contains_source_type(data, "temporary_demo"):
        raise ConfigurationError("scientific mode forbids temporary_demo parameters.")
    if len(data["species"]) != 15 or set(data["species"]) != {f"C{i}" for i in range(1, 16)}:
        raise ConfigurationError("Exactly C1...C15 species are required.")
    reactions = data["kinetics"]["reactions"]
    if len(reactions) != 15 or set(reactions) != {f"R{i}" for i in range(1, 16)}:
        raise ConfigurationError("Exactly R1...R15 reactions are required.")
    if reactions["R15"].get("reversible") is not False or "reverse" in reactions["R15"]:
        raise ConfigurationError("R15 must be irreversible and must not contain a reverse block.")
    for i in range(1, 16):
        reaction = reactions[f"R{i}"]
        if "forward" not in reaction:
            raise ConfigurationError(f"R{i} is missing forward parameters.")
        if reaction.get("reversible") and "reverse" not in reaction:
            raise ConfigurationError(f"R{i} is reversible but has no reverse parameters.")
    _positive(data["conditions"]["temperature"]["value"], "temperature")
    _positive(data["conditions"]["pressure"]["value"], "pressure")
    _positive(data["pellet"]["radius"]["value"], "pellet.radius")
    porosity = float(data["pellet"]["porosity"]["value"])
    if not 0.0 < porosity < 1.0:
        raise ConfigurationError("pellet.porosity must be between 0 and 1.")
    _positive(data["pellet"]["tortuosity"]["value"], "pellet.tortuosity")
    _positive(data["pellet"]["pore_radius"]["value"], "pellet.pore_radius")
    nr = int(data["numerics"]["spatial"]["radial_cells"])
    if nr < 2:
        raise ConfigurationError("radial_cells must be at least 2.")
    if data["diffusion"]["model"] != "knudsen":
        raise ConfigurationError("The first V1 implementation currently supports diffusion.model: knudsen.")
    initial_model = data["initial_conditions"]["model"]
    if initial_model not in {"bulk_equilibrated", "hydrogen_pretreated", "user_defined"}:
        raise ConfigurationError("Unsupported initial_conditions.model.")
    if initial_model == "user_defined" and "concentrations" not in data["initial_conditions"]:
        raise ConfigurationError("user_defined initial state requires a concentrations mapping.")
    simulation_mode = data["numerics"]["time"]["mode"]
    if simulation_mode not in {"transient", "steady"}:
        raise ConfigurationError("numerics.time.mode must be transient or steady.")
    if data["numerics"]["solver"]["linear_solver"] != "dense":
        raise ConfigurationError("The current CVODE backend supports linear_solver: dense.")
    for key in ("C1", "C4"):
        _nonnegative(data["bulk"]["concentrations"][key]["value"], f"bulk.{key}")
    if float(data["bulk"]["concentrations"]["C1"]["value"]) == 0.0:
        raise ConfigurationError("C1 bulk concentration must be positive to calculate conversion.")
    if data["bulk"].get("composition_input") == "mole_fraction":
        fractions = [float(data["bulk"]["mole_fractions"][f"C{i}"]) for i in range(1, 16)]
        if any(value < 0.0 for value in fractions):
            raise ConfigurationError("Bulk mole fractions must be nonnegative.")
        if not math.isclose(sum(fractions), 1.0, rel_tol=0.0, abs_tol=1.0e-10):
            raise ConfigurationError("Bulk mole fractions must sum to one.")


def normalized_hash(data: dict[str, Any]) -> str:
    payload = json.dumps(data, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]


def dump_normalized(data: dict[str, Any], path: Path) -> None:
    with path.open("w", encoding="utf-8") as stream:
        yaml.safe_dump(data, stream, allow_unicode=True, sort_keys=False)


def write_fortran_namelist(data: dict[str, Any], path: Path) -> None:
    species = [data["species"][f"C{i}"] for i in range(1, 16)]
    reactions = [data["kinetics"]["reactions"][f"R{i}"] for i in range(1, 16)]
    bulk = [float(data["bulk"]["concentrations"][f"C{i}"]["value"]) for i in range(1, 16)]
    initial = _initial_concentrations(data, bulk)
    molar_mass = [float(item["molar_mass"]["value"]) for item in species]
    ln_f = [float(item["forward"]["ln_k0"]["value"]) for item in reactions]
    e_f = [_energy_j_per_mol(item["forward"]["E"]) for item in reactions]
    ln_r = [float(item.get("reverse", {}).get("ln_k0", {}).get("value", 0.0)) for item in reactions]
    e_r = [_energy_j_per_mol(item["reverse"]["E"]) if "reverse" in item else 0.0 for item in reactions]
    reversible = [bool(item["reversible"]) for item in reactions]
    numerics = data["numerics"]
    solver = numerics["solver"]
    time = numerics["time"]
    lines = ["&simulation"]
    scalar_values = {
        "run_mode": _quote(data["run_mode"]),
        "run_name": _quote(data["run_name"]),
        "simulation_mode": _quote(time["mode"]),
        "nr": int(numerics["spatial"]["radial_cells"]),
        "temperature": float(data["conditions"]["temperature"]["value"]),
        "pressure": float(data["conditions"]["pressure"]["value"]),
        "radius": float(data["pellet"]["radius"]["value"]),
        "porosity": float(data["pellet"]["porosity"]["value"]),
        "tortuosity": float(data["pellet"]["tortuosity"]["value"]),
        "pore_radius": float(data["pellet"]["pore_radius"]["value"]),
        "acid_activity": float(data["kinetics"]["acid_sites"]["activity"]["value"]),
        "concentration_reference": float(data["kinetics"]["demo_embedding"]["normalization_concentration"]["value"]),
        "t_start": float(time["start_time"]["value"]),
        "t_end": float(time["end_time"]["value"]),
        "output_interval": float(time["output_interval"]["value"]),
        "rtol": float(solver["rtol"]),
        "atol": float(solver["atol"]["value"]),
        "initial_step": float(solver["initial_step"]["value"]),
        "max_step": float(solver["max_step"]["value"]),
        "max_steps": int(solver["max_steps"]),
        "negative_warning": float(data["nonnegativity"]["warning_tolerance"]["value"]),
        "negative_error": float(data["nonnegativity"]["error_tolerance"]["value"]),
        "steady_tolerance": float(numerics["steady_state"]["derivative_tolerance"]["value"]),
        "steady_maximum_time": float(numerics["steady_state"]["maximum_time"]["value"]),
    }
    for key, value in scalar_values.items():
        lines.append(f"  {key} = {_fortran_scalar(value)}")
    for key, values in {
        "molar_mass": molar_mass,
        "c_bulk": bulk,
        "c_initial": initial,
        "ln_k0_forward": ln_f,
        "activation_forward": e_f,
        "ln_k0_reverse": ln_r,
        "activation_reverse": e_r,
        "reversible": reversible,
    }.items():
        lines.append(f"  {key} = {', '.join(_fortran_scalar(value) for value in values)}")
    lines.append("/")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _energy_j_per_mol(parameter: dict[str, Any]) -> float:
    value = float(parameter["value"])
    unit = parameter["unit"]
    if unit == "J/mol":
        return value
    if unit == "kJ/mol":
        return value * 1000.0
    if unit == "kcal/mol":
        return value * 4184.0
    raise ConfigurationError(f"Unsupported activation-energy unit: {unit}")


def _initial_concentrations(data: dict[str, Any], bulk: list[float]) -> list[float]:
    model = data["initial_conditions"]["model"]
    if model == "bulk_equilibrated":
        return bulk.copy()
    if model == "hydrogen_pretreated":
        values = [0.0]*15
        values[3] = bulk[3]
        return values
    mapping = data["initial_conditions"]["concentrations"]
    return [float(mapping[f"C{i}"]["value"]) for i in range(1, 16)]


def _contains_source_type(node: Any, source_type: str) -> bool:
    if isinstance(node, dict):
        if node.get("source_type") == source_type:
            return True
        return any(_contains_source_type(value, source_type) for value in node.values())
    if isinstance(node, list):
        return any(_contains_source_type(value, source_type) for value in node)
    return False


def _positive(value: Any, name: str) -> None:
    if not isinstance(value, (int, float)) or isinstance(value, bool) or value <= 0:
        raise ConfigurationError(f"{name} must be a positive number.")


def _nonnegative(value: Any, name: str) -> None:
    if not isinstance(value, (int, float)) or isinstance(value, bool) or value < 0:
        raise ConfigurationError(f"{name} must be a nonnegative number.")


def _quote(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def _fortran_scalar(value: Any) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:.17e}".replace("e", "d")
    raise TypeError(f"Cannot serialize {value!r} to Fortran namelist")

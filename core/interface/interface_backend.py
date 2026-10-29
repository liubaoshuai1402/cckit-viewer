"""Engine dispatch and validated settings for the Interface Builder feature."""

from __future__ import annotations

import copy
import math
from pathlib import Path
from typing import Any, Callable, Mapping


SEARCH_KEYS = (
    "engine", "film_structure", "substrate_plane", "substrate_plane_notation",
    "film_plane", "film_plane_notation", "max_area", "max_area_ratio_tol",
    "max_length_tol", "max_angle_tol", "bidirectional", "lowest", "max_matches",
    "film_termination_ftol", "substrate_termination_ftol", "label_index",
    "filter_out_sym_slabs",
    "im_angle_start", "im_angle_stop", "im_angle_step", "im_max_denominator",
    "im_max_strain", "im_max_area_change", "im_max_sigma", "im_max_index",
    "im_rank_by", "im_rotation_tol",
)


def default_config() -> dict[str, Any]:
    """Use pymatgen's ZSL defaults, with explicit termination enumeration."""

    return {
        "mode": "build", "engine": "pymatgen", "film_structure": "", "system": "",
        "substrate_plane": [0, 0, 1], "substrate_plane_notation": "three",
        "film_plane": [0, 0, 1], "film_plane_notation": "three",
        "max_area": 400.0, "max_area_ratio_tol": 0.09,
        "max_length_tol": 0.03, "max_angle_tol": 0.01,
        "bidirectional": False, "lowest": False, "max_matches": 100,
        "film_termination_ftol": 0.25, "substrate_termination_ftol": 0.25,
        "label_index": True, "filter_out_sym_slabs": True,
        "match_index": 0, "termination": None,
        "film_thickness": 3.0, "substrate_thickness": 3.0, "in_layers": True,
        "gap": 2.0, "vacuum_over_film": 20.0,
        "shift_mode": "custom", "custom_shift_basis": "supercell",
        "offset_a": 0.0, "offset_b": 0.0, "offset_grid": [1, 1],
        "im_angle_start": 0.0, "im_angle_stop": 60.0, "im_angle_step": 1.0,
        "im_max_denominator": 30, "im_max_strain": 0.03,
        "im_max_area_change": 0.05, "im_max_sigma": 100.0,
        "im_max_index": 400, "im_rank_by": "sigma",
        "im_rotation_tol": 0.5,
        "im_shift_basis": "dscl", "im_contact_depth": 1.0,
        "im_registry_tolerance": 0.3, "im_min_distance": 1.0,
        "im_max_atoms": 50000,
    }


def parameters(config: Mapping[str, Any]) -> dict[str, Any]:
    """Validate controls independently of Qt and optional modeling imports."""

    if not isinstance(config, Mapping):
        raise ValueError("interface.json must contain a JSON object")
    if "engine" not in config and any(key in config for key in ("csl_mode", "cnid_grid")):
        raise ValueError("Legacy InterfaceMaster JSON: create a new configuration with an explicit engine and current settings.")
    values = default_config()
    values.update(copy.deepcopy(dict(config)))
    if values["engine"] not in {"pymatgen", "interface-master"}:
        raise ValueError("engine must be pymatgen or interface-master")
    if values["mode"] != "build":
        raise ValueError("Interface Builder supports build mode only")
    system = values["system"]
    if not isinstance(system, str) or len(system) > 200 or any(not char.isprintable() for char in system):
        raise ValueError("system must be a single-line label of at most 200 printable characters")
    values["system"] = system.strip()
    for key in ("film_structure", "substrate_structure"):
        if key in values and not isinstance(values[key], str):
            raise ValueError(f"{key} must be a file path")
    for side in ("substrate", "film"):
        key = f"{side}_plane"
        indices = values[key]
        if not isinstance(indices, (list, tuple)) or len(indices) not in {3, 4}:
            raise ValueError(f"{key} needs three or four integer indices")
        if any(isinstance(i, bool) or not isinstance(i, int) for i in indices):
            raise ValueError(f"{key} must contain integer indices")
        notation = config.get(f"{key}_notation", "four" if len(indices) == 4 else "three")
        if notation not in {"three", "four"} or len(indices) != (4 if notation == "four" else 3):
            raise ValueError(f"{key}: index count does not match notation")
        if len(indices) == 4 and sum(indices[:3]) != 0:
            raise ValueError(f"{key}: hexagonal indices require h + k + i = 0")
        if not any(indices):
            raise ValueError(f"{key} cannot be zero")
        values[key] = list(indices)
        values[f"{key}_notation"] = notation
    for key in ("bidirectional", "lowest", "label_index", "filter_out_sym_slabs", "in_layers"):
        if not isinstance(values[key], bool):
            raise ValueError(f"{key} must be a JSON boolean")
    for key in ("max_matches", "match_index", "im_max_denominator", "im_max_index", "im_max_atoms"):
        value = values[key]
        minimum = 0 if key == "match_index" else 1
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            raise ValueError(f"{key} must be an integer >= {minimum}")
        maximum = {"max_matches": 10000, "match_index": 9999,
                   "im_max_denominator": 100, "im_max_index": 10000,
                   "im_max_atoms": 200000}[key]
        if value > maximum:
            raise ValueError(f"{key} exceeds the supported candidate range")
    positive = {"max_area", "film_thickness", "substrate_thickness",
                "film_termination_ftol", "substrate_termination_ftol"}
    for key in (
        "max_area", "max_area_ratio_tol", "max_length_tol", "max_angle_tol",
        "film_termination_ftol", "substrate_termination_ftol", "film_thickness",
        "substrate_thickness", "gap", "vacuum_over_film", "offset_a", "offset_b",
        "im_angle_start", "im_angle_stop", "im_angle_step", "im_max_strain",
        "im_max_area_change", "im_max_sigma", "im_contact_depth",
        "im_registry_tolerance", "im_min_distance",
        "im_rotation_tol",
    ):
        raw = values[key]
        if isinstance(raw, bool):
            raise ValueError(f"{key} must be a finite number")
        try:
            value = float(raw)
        except (ValueError, TypeError) as exc:
            raise ValueError(f"{key} must be a finite number") from exc
        if not math.isfinite(value):
            raise ValueError(f"{key} must be a finite number")
        if key not in {"offset_a", "offset_b", "im_angle_start", "im_angle_stop"} and (value < 0 or (key in positive and value == 0)):
            raise ValueError(f"{key} is outside its allowed range")
        values[key] = value
    # Keep JSON and the visible control ranges consistent instead of silently
    # clamping loaded values while building from different hidden settings.
    bounds = {
        "max_area": (1, 100000),
        "max_area_ratio_tol": (0, 1), "max_length_tol": (0, 1), "max_angle_tol": (0, 1),
        "film_termination_ftol": (0.000001, 5), "substrate_termination_ftol": (0.000001, 5),
        "film_thickness": (1 if values["in_layers"] else 0.0001, 10000),
        "substrate_thickness": (1 if values["in_layers"] else 0.0001, 10000),
        "gap": (0, 10000), "vacuum_over_film": (0, 10000),
        "offset_a": (-100, 100), "offset_b": (-100, 100),
        "im_angle_start": (-360, 360), "im_angle_stop": (-360, 360),
        "im_angle_step": (0.001, 360), "im_max_strain": (0, 0.5),
        "im_max_area_change": (0, 1), "im_max_sigma": (1, 10000),
        "im_contact_depth": (0, 100), "im_registry_tolerance": (0, 10),
        "im_min_distance": (0, 10),
        "im_rotation_tol": (0, 10),
    }
    for key, (lower, upper) in bounds.items():
        if not lower <= values[key] <= upper:
            raise ValueError(f"{key} must be between {lower} and {upper}")
    if values["im_angle_stop"] < values["im_angle_start"]:
        raise ValueError("Rotation stop must be greater than or equal to start")
    angle_count = math.ceil((values["im_angle_stop"] - values["im_angle_start"]) / values["im_angle_step"]) + 1
    if angle_count * values["im_max_denominator"] > 250000:
        raise ValueError("Rotation/denominator search exceeds 250000 trials; increase the angle step or reduce the bounds")
    if values["im_rank_by"] not in {"sigma", "strain", "area"}:
        raise ValueError("im_rank_by must be sigma, strain or area")
    if values["im_shift_basis"] not in {"dscl", "supercell"}:
        raise ValueError("im_shift_basis must be dscl or supercell")
    if values["custom_shift_basis"] not in {"dscl", "supercell"}:
        raise ValueError("custom_shift_basis must be dscl or supercell")
    if values["in_layers"]:
        for key in ("film_thickness", "substrate_thickness"):
            if not values[key].is_integer():
                raise ValueError("Thickness in unit planes must be a positive integer")
    grid = values["offset_grid"]
    if not isinstance(grid, (list, tuple)) or len(grid) != 2 or any(
        isinstance(n, bool) or not isinstance(n, int) or not 1 <= n <= 100 for n in grid
    ):
        raise ValueError("offset_grid requires two integers between 1 and 100")
    values["offset_grid"] = list(grid)
    if "shift_mode" not in config:
        sampling = values["offset_grid"] != [1, 1]
        shifted = any(values[key] % 1 != 0 for key in ("offset_a", "offset_b"))
        if sampling and shifted:
            raise ValueError(
                "This older configuration combines a custom shift and sampling. "
                "Set shift_mode to 'custom' or 'sampling' explicitly."
            )
        values["shift_mode"] = "sampling" if sampling else "custom"
    if values["shift_mode"] not in {"custom", "sampling"}:
        raise ValueError("shift_mode must be custom or sampling")
    termination = values["termination"]
    if termination is not None and (
        not isinstance(termination, (list, tuple)) or len(termination) != 2
        or any(not isinstance(label, str) for label in termination)
    ):
        raise ValueError("termination must be null or [film label, substrate label]")
    if termination is not None:
        values["termination"] = list(termination)
    return values


def shift_grid(config: Mapping[str, Any]) -> tuple[int, int]:
    """Return the active grid; saved sampling settings do not affect custom shifts."""

    if config["shift_mode"] == "custom":
        return (1, 1)
    n1, n2 = config["offset_grid"]
    return n1, n2


def frame_count(config: Mapping[str, Any]) -> int:
    """Use the same active frame count for preview construction and export."""

    return math.prod(shift_grid(config))


def active_shift_basis(config: Mapping[str, Any]) -> str:
    """Keep custom and sampling bases independent, with a/b only for ZSL."""

    if config["engine"] != "interface-master":
        return "supercell"
    if config["shift_mode"] == "custom":
        return config.get("custom_shift_basis", "supercell")
    return config["im_shift_basis"]


def frame_identity(
    config: Mapping[str, Any], system: str, coefficients: list[float], record: Mapping[str, Any],
) -> dict[str, Any]:
    """Format the XYZ tag while retaining unrounded frame values for calculations."""

    def number(value: float) -> str:
        value = float(value)
        return format(0.0 if value == 0 else value, ".12g")

    system = config["system"] or system
    d1, d2 = (float(value) for value in coefficients)
    vacuum = float(config["vacuum_over_film"])
    zspacing = float(config["gap"])
    if config["engine"] == "interface-master":
        matching = f"sigma_{record['sigma']:.2f}_theta_{record['twist_angle_deg']:.1f}"
    else:
        matching = f"zsl_match_{record['index'] + 1}_strain_{number(record['von_mises_strain'])}"
    return {
        "system": system,
        "config_type": f"{system}_{matching}_zspacing_{zspacing:.1f}_d1_{d1:.1f}_d2_{d2:.1f}",
        "d1": d1, "d2": d2, "shift_basis": active_shift_basis(config),
        "vacuum_A": vacuum, "gap_A": zspacing, "zspacing_A": zspacing,
    }


def _engine(config: Mapping[str, Any]) -> Any:
    if config["engine"] == "interface-master":
        from core.interface import interface_master_engine

        return interface_master_engine
    from core.interface import pymatgen_engine

    return pymatgen_engine


def search_matches(
    input_path: Path, config: Mapping[str, Any], *, cancelled: Callable[[], bool] | None = None,
) -> dict[str, Any]:
    """Return matching cells and termination choices for the active engine."""

    values = parameters(config)
    return _engine(values).search_matches(input_path, values, cancelled=cancelled)


def build_interface_frames(
    input_path: Path, config: Mapping[str, Any], *,
    progress: Callable[[int, int], None] | None = None,
    cancelled: Callable[[], bool] | None = None,
) -> dict[str, Any]:
    """Build the selected match using the active engine."""

    values = parameters(config)
    return _engine(values).build_frames(input_path, values, progress=progress, cancelled=cancelled)


def _default_output_name(mode: str, config: Mapping[str, Any]) -> str:
    """Give exported structures a feature-specific basename."""

    return "POSCAR_Interface"

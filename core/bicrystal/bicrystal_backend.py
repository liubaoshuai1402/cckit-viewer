# CCkit 功能说明：
# 扫描旋转轴的 CSL 候选，并按 gbin.json 构建对称倾斜晶界 POSCAR。
"""Bicrystal Builder: build rotational symmetric-tilt boundaries with InterfaceMaster.

The command has two modes selected by ``gbin.json`` (or ``--mode``):

``scan``
    Sweep a rotation-angle range around a crystallographic axis and write a
    CSV containing only angles for which a CSL candidate was found within the
    configured maximum sigma.

``build``
    Select one angle/sigma pair, convert the requested symmetric GB plane to
    the primitive basis used by InterfaceMaster, and write a VASP POSCAR.

The CSL policy is configured with ``csl_mode``: ``"exact"`` calls
``search_fixed(..., exact=True)`` and forces ``D=I`` for an exact CSL.  In
that mode InterfaceMaster may adjust the rotation to its exact CSL value.
``"near"`` (the default) keeps the requested rotation and permits the small
deformation matrix used by InterfaceMaster for a near-CSL construction.  The
legacy ``exact`` boolean key remains supported.

Axis and plane indices refer to the basis of the input structure.
Standard hexagonal four-index notation is accepted through
``rotation_axis_notation`` and ``gb_plane_notation`` ("three" or "four").
Interface translations are Cartesian Angstrom vectors; ``dx`` is the normal
spacing in Angstrom, while ``dp1`` and ``dp2`` select termination shifts.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import json
import math
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence


def default_viewer_config() -> dict[str, Any]:
    """Start the GUI in scan mode with independent Cu-example parameters."""

    return {
        "mode": "scan",
        "csl_mode": "near",
        "rotation_axis": [0, 0, 1],
        "rotation_angle_deg": 32.0,
        "sigma": 5,
        "max_sigma": 6,
        "gb_plane": [2, -1, 0],
        "mirror": True,
        "normal_distance": 0.0,
        "vacuum": 0.0,
        "normal_ortho": False,
        "plane_ortho": False,
        "termination_shift": {"dp1": 0.0, "dp2": 0.0},
        "tolerances": {
            "du": 0.2,
            "area": 0.2,
            "deformation": 0.1,
            "integer": 0.01,
        },
        "angle_start_deg": 0.0,
        "angle_stop_deg": 180.0,
        "angle_step_deg": 1.0,
        "cnid_grid": [4, 4],
    }


def _backend() -> tuple[Any, Any, Any, Any, Any, Any]:
    """Import optional InterfaceMaster dependencies with an actionable error."""

    try:
        import numpy as np
        from pymatgen.core import Structure
        from interfacemaster.cellcalc import get_primitive_hkl, rot
        from interfacemaster.interface_generator import (
            convert_vector_index,
            core,
        )
    except ImportError as exc:
        raise RuntimeError(
            "model-gb requires numpy, pymatgen, and interfacemaster; "
            "install InterfaceMaster in the Python environment used by CCkit"
        ) from exc

    return np, Structure, get_primitive_hkl, rot, convert_vector_index, core


def index_notation(value: Any, notation: str | None = None) -> str:
    """Choose three/four indices, inferring notation for legacy JSON files."""

    if notation is None:
        return "four" if isinstance(value, (list, tuple)) and len(value) == 4 else "three"
    if notation not in {"three", "four"}:
        raise ValueError("index notation must be 'three' or 'four'")
    return notation


def index_values(value: Any, kind: str, notation: str | None = None) -> list[float]:
    """Validate direction [uvw]/[uvtw] or plane (hkl)/(hkil) indices."""

    notation = index_notation(value, notation)
    length = 4 if notation == "four" else 3
    if not isinstance(value, (list, tuple)) or len(value) != length:
        raise ValueError(f"{kind} requires {length} indices in {notation}-index notation")
    if any(isinstance(item, bool) for item in value):
        raise ValueError(f"{kind} indices must be numbers")
    values = [_finite_float(item, f"{kind} index") for item in value]
    if not any(values):
        raise ValueError(f"{kind} indices must not all be zero")
    if kind not in {"direction", "plane"}:
        raise ValueError("index kind must be 'direction' or 'plane'")
    if notation == "four":
        if any(item != round(item) for item in values):
            raise ValueError("four-index notation requires integer indices")
        if abs(sum(values[:3])) > 1e-10:
            constraint = "u+v+t=0" if kind == "direction" else "h+k+i=0"
            raise ValueError(f"four-index {kind} must satisfy {constraint}")
    return values


def three_indices(
    value: Any, kind: str, notation: str | None = None, lattice: Any = None
) -> list[float]:
    """Convert standard hexagonal four-index notation into a three-axis basis.

    Lattice vectors, when provided, are Cartesian columns. Four indices require
    a1/a2 of equal length at 120 degrees and c normal to the basal plane. The
    function never silently standardizes or changes the supplied unit cell.
    """

    notation = index_notation(value, notation)
    values = index_values(value, kind, notation)
    if notation == "three":
        return values
    if lattice is not None:
        import numpy as np

        basis = np.asarray(lattice, dtype=float)
        lengths = np.linalg.norm(basis, axis=0)
        if basis.shape != (3, 3) or not np.all(np.isfinite(basis)) or np.any(lengths <= 0):
            raise ValueError("four-index notation requires a valid hexagonal lattice")
        unit = basis / lengths
        if not (
            abs(lengths[0] / lengths[1] - 1) < 1e-3
            and abs(np.dot(unit[:, 0], unit[:, 1]) + 0.5) < 2e-3
            and abs(np.dot(unit[:, 0], unit[:, 2])) < 2e-3
            and abs(np.dot(unit[:, 1], unit[:, 2])) < 2e-3
        ):
            raise ValueError(
                "four-index notation requires a standard hexagonal input cell: "
                "a=b, gamma=120 degrees, c perpendicular to a/b. "
                "Use three indices for the current basis or supply a standard hexagonal cell."
            )
    first, second, redundant, last = values
    if kind == "plane":
        return [first, second, last]
    converted = [int(first - redundant), int(second - redundant), int(last)]
    divisor = math.gcd(*converted)
    return [item // divisor for item in converted]


def four_indices(value: Any, kind: str) -> list[int]:
    """Express integral three-axis indices in standard hexagonal notation."""

    first, second, last = index_values(value, kind, "three")
    if any(item != round(item) for item in (first, second, last)):
        raise ValueError("switching to four indices requires integral three-index input")
    first, second, last = int(first), int(second), int(last)
    if kind == "plane":
        return [first, second, -first - second, last]
    converted = [2 * first - second, 2 * second - first, -first - second, 3 * last]
    divisor = math.gcd(*converted)
    return [item // divisor for item in converted]


def _finite_vector(
    value: Any,
    name: str,
    length: int,
    np: Any,
    *,
    allow_two: bool = False,
) -> Any:
    """Validate and convert a numeric vector from the JSON configuration."""

    if not isinstance(value, (list, tuple)):
        raise ValueError(f"{name} must be a list of {length} numbers")
    if allow_two and len(value) == 2:
        value = [0.0, value[0], value[1]]
    if len(value) != length:
        raise ValueError(f"{name} must contain {length} numbers")
    try:
        vector = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain only numbers") from exc
    if not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must contain only finite numbers")
    return vector


def _positive_int(value: Any, name: str) -> int:
    """Return a strictly positive integer configuration value."""

    if isinstance(value, bool):
        raise ValueError(f"{name} must be a positive integer")
    try:
        number = int(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a positive integer") from exc
    if number <= 0 or float(value) != number:
        raise ValueError(f"{name} must be a positive integer")
    return number


def _nonnegative_float(value: Any, name: str) -> float:
    """Return a finite non-negative float configuration value."""

    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a non-negative number") from exc
    if not math.isfinite(number) or number < 0:
        raise ValueError(f"{name} must be a non-negative number")
    return number


def _finite_float(value: Any, name: str) -> float:
    """Return a finite float, allowing signed termination shifts."""

    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a number") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _repeat_vector(value: Any, name: str) -> list[int]:
    """Validate a three-component positive supercell repeat vector."""

    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError(f"{name} must contain three positive integers")
    return [_positive_int(component, name) for component in value]


def _config_value(
    config: Mapping[str, Any],
    names: Sequence[str],
    default: Any = None,
) -> Any:
    """Return the first present value from a list of compatible key names."""

    for name in names:
        if name in config:
            return config[name]
    return default


def _read_json(path: Path) -> dict[str, Any]:
    """Read a JSON object and report malformed configuration clearly."""

    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid JSON in {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise ValueError(f"configuration root must be a JSON object: {path}")
    return data


def _resolve_input(path_value: str) -> Path:
    """Resolve and validate the structure input path."""

    path = Path(path_value).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"input structure was not found: {path}")
    return path


def _resolve_config(path_value: str | None, input_path: Path) -> Path:
    """Resolve gbin.json, defaulting to the input structure directory."""

    if path_value is None:
        candidates = (input_path.parent / "gbin.json", Path.cwd() / "gbin.json")
        path = next(
            (candidate for candidate in candidates if candidate.is_file()),
            candidates[0],
        )
    else:
        path = Path(path_value).expanduser()
        if not path.is_absolute():
            path = Path.cwd() / path
    path = path.resolve()
    if not path.is_file():
        if path_value is None:
            raise FileNotFoundError(
                "configuration file was not found; expected gbin.json beside "
                f"the input ({candidates[0]}) or in the current directory "
                f"({candidates[1]})"
            )
        raise FileNotFoundError(f"configuration file was not found: {path}")
    return path


def _resolve_output(
    path_value: str | None,
    input_path: Path,
    default_name: str,
) -> Path:
    """Apply CCkit's input-directory rule to a structure output path."""

    if path_value is None:
        path = input_path.parent / default_name
    else:
        path = Path(path_value).expanduser()
        if not path.is_absolute():
            if len(path.parts) == 1:
                path = input_path.parent / path
            else:
                path = Path.cwd() / path
    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _ensure_poscar_prefix(path: Path) -> Path:
    """Ensure generated bicrystals use a POSCAR_ filename prefix."""

    if path.name.startswith("POSCAR_"):
        return path
    return path.with_name(f"POSCAR_{path.name}")


def _tolerances(config: Mapping[str, Any]) -> dict[str, float]:
    """Read InterfaceMaster matching tolerances from nested or flat keys."""

    nested = config.get("tolerances", {})
    if not isinstance(nested, dict):
        raise ValueError("tolerances must be a JSON object")

    values = {
        "du": _config_value(config, ("du",), nested.get("du", 0.02)),
        "area": _config_value(
            config,
            ("area_tolerance", "area_tol"),
            nested.get("area", 0.02),
        ),
        "deformation": _config_value(
            config,
            ("deformation_tolerance", "deformation_tol"),
            nested.get("deformation", 0.02),
        ),
        "integer": _config_value(
            config,
            ("integer_tolerance", "integer_tol"),
            nested.get("integer", 1.0e-3),
        ),
    }
    result = {}
    for name, value in values.items():
        result[name] = _nonnegative_float(value, f"tolerances.{name}")
        if result[name] == 0:
            raise ValueError(f"tolerances.{name} must be greater than zero")
    return result


def _csl_mode(config: Mapping[str, Any]) -> tuple[str, bool]:
    """Read and validate exact-CSL versus near-CSL selection."""

    mode_value = config.get("csl_mode")
    legacy_value = config.get("exact")
    if mode_value is not None:
        if not isinstance(mode_value, str):
            raise ValueError("csl_mode must be either 'exact' or 'near'")
        normalized = mode_value.strip().lower().replace("_", "-")
        aliases = {
            "exact": ("exact", True),
            "near": ("near", False),
            "near-csl": ("near", False),
            "approximate": ("near", False),
        }
        if normalized not in aliases:
            raise ValueError("csl_mode must be either 'exact' or 'near'")
        mode, exact = aliases[normalized]
        if legacy_value is not None and _parse_bool(legacy_value, "exact") != exact:
            raise ValueError("csl_mode and exact specify conflicting CSL policies")
        return mode, exact

    if legacy_value is None:
        return "near", False
    exact = _parse_bool(legacy_value, "exact")
    return ("exact", True) if exact else ("near", False)


def _parse_bool(value: Any, name: str) -> bool:
    """Parse JSON booleans and common textual boolean spellings."""

    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "yes", "on", "1"}:
            return True
        if normalized in {"false", "no", "off", "0"}:
            return False
    raise ValueError(f"{name} must be a boolean")


def _axis_and_rotation(
    structure: Any,
    axis_input: Any,
    angle_deg: float,
    np: Any,
    rot: Any,
    convert_vector_index: Any,
    core: Any,
) -> tuple[Any, Any, Any, Any]:
    """Convert an input-basis axis and construct the requested rotation."""

    axis_conventional = _finite_vector(axis_input, "rotation_axis", 3, np)
    if np.linalg.norm(axis_conventional) == 0:
        raise ValueError("rotation_axis must not be the zero vector")

    base = core(structure, structure, verbose=False)
    axis_primitive = convert_vector_index(
        base.conv_lattice_1,
        base.lattice_1,
        axis_conventional,
    )
    axis_cartesian = np.dot(base.lattice_1, axis_primitive)
    axis_cartesian = axis_cartesian / np.linalg.norm(axis_cartesian)
    requested_rotation = rot(axis_cartesian, np.deg2rad(angle_deg))
    return axis_conventional, axis_primitive, axis_cartesian, requested_rotation


def _new_interface(
    structure: Any,
    max_sigma: int,
    tolerances: Mapping[str, float],
    core: Any,
) -> Any:
    """Create and configure an InterfaceMaster core object."""

    interface = core(structure, structure, verbose=False)
    interface.parse_limit(
        tolerances["du"],
        tolerances["area"],
        max_sigma,
        max_sigma,
        tolerances["deformation"],
    )
    # search_fixed() does not return a status.  None lets us distinguish a
    # successful search from its default identity matrix after the call.
    interface.R = None
    return interface


def _search_angle(
    structure: Any,
    axis_input: Any,
    angle_deg: float,
    max_sigma: int,
    tolerances: Mapping[str, float],
    exact: bool,
    np: Any,
    rot: Any,
    convert_vector_index: Any,
    core: Any,
) -> tuple[Any, dict[str, Any]]:
    """Search one angle and return its InterfaceMaster object and CSV row."""

    _, _, axis_cartesian, requested_rotation = _axis_and_rotation(
        structure,
        axis_input,
        angle_deg,
        np,
        rot,
        convert_vector_index,
        core,
    )
    interface = _new_interface(structure, max_sigma, tolerances, core)

    # InterfaceMaster writes log.fixed_search and several build intermediates
    # relative to the process directory.  Keep those files in a disposable
    # directory so CCkit never leaves calculation artifacts beside the input.
    with tempfile.TemporaryDirectory(prefix="cckit_gb_") as work_dir:
        old_cwd = Path.cwd()
        try:
            os.chdir(work_dir)
            with contextlib.redirect_stdout(io.StringIO()):
                interface.search_fixed(
                    requested_rotation,
                    exact=exact,
                    tol=tolerances["integer"],
                )
        finally:
            os.chdir(old_cwd)

    row: dict[str, Any] = {
        "angle_deg": angle_deg,
        "sigma_1": "",
        "sigma_2": "",
        "deformation_max_abs": "",
        "status": "not_found",
    }
    if interface.R is None:
        return interface, row

    sigma_1 = int(round(abs(np.linalg.det(interface.U1))))
    sigma_2 = int(round(abs(np.linalg.det(interface.U2))))
    deformation_max_abs = float(
        np.max(np.abs(np.asarray(interface.D, dtype=float) - np.eye(3)))
    )
    row.update(
        sigma_1=sigma_1,
        sigma_2=sigma_2,
        deformation_max_abs=deformation_max_abs,
        status="found",
    )
    return interface, row


def _scan_angles(start: float, stop: float, step: float) -> list[float]:
    """Generate an inclusive, numerically stable angle grid."""

    if step <= 0:
        raise ValueError("angle_step_deg must be greater than zero")
    if stop < start:
        raise ValueError(
            "angle_stop_deg must be greater than or equal to angle_start_deg"
        )
    count = int(math.floor((stop - start) / step + 1.0e-10)) + 1
    return [start + index * step for index in range(count)]


def scan_candidates(
    input_path: Path,
    config: Mapping[str, Any],
    output_path: Path,
) -> int:
    """Scan a rotation axis and write all tested CSL candidates to CSV."""

    np, Structure, _, rot, convert_vector_index, core = _backend()
    structure = Structure.from_file(input_path)
    axis = _config_value(config, ("rotation_axis", "axis"), None)
    if axis is None:
        raise ValueError("scan mode requires rotation_axis, for example [0, 0, 1]")
    axis = three_indices(
        axis, "direction", config.get("rotation_axis_notation"), structure.lattice.matrix.T
    )
    max_sigma = _positive_int(
        _config_value(config, ("max_sigma", "maxsigma"), None),
        "max_sigma",
    )
    tolerances = _tolerances(config)
    csl_mode, exact = _csl_mode(config)
    angle_range = _config_value(
        config,
        ("angle_range_deg", "angle_range", "theta_range"),
        None,
    )
    if angle_range is not None:
        if not isinstance(angle_range, (list, tuple)) or len(angle_range) != 2:
            raise ValueError("angle_range must contain [start_deg, stop_deg]")
        range_start, range_stop = angle_range
    else:
        range_start, range_stop = 0.0, 180.0
    start = float(config.get("angle_start_deg", range_start))
    stop = float(config.get("angle_stop_deg", range_stop))
    step = float(config.get("angle_step_deg", config.get("angle_step", 1.0)))
    for value, name in (
        (start, "angle_start_deg"),
        (stop, "angle_stop_deg"),
        (step, "angle_step_deg"),
    ):
        if not math.isfinite(value):
            raise ValueError(f"{name} must be finite")

    rows = []
    for angle_deg in _scan_angles(start, stop, step):
        _, row = _search_angle(
            structure,
            axis,
            angle_deg,
            max_sigma,
            tolerances,
            exact,
            np,
            rot,
            convert_vector_index,
            core,
        )
        row["csl_mode"] = csl_mode
        rows.append(row)

    found_rows = [row for row in rows if row["status"] == "found"]

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "angle_deg",
                "sigma_1",
                "sigma_2",
                "deformation_max_abs",
                "csl_mode",
                "status",
            ),
        )
        writer.writeheader()
        writer.writerows(found_rows)

    found = len(found_rows)
    print(
        f"Scanned {len(rows)} angles in {csl_mode}-CSL mode; "
        f"found {found} CSL candidates."
    )
    print(f"CSV: {output_path}")
    return found


def _termination_values(config: Mapping[str, Any]) -> tuple[float, float]:
    """Read the two termination shifts used by InterfaceMaster."""

    shifts = config.get("termination_shift", {})
    if shifts is None:
        shifts = {}
    if not isinstance(shifts, dict):
        raise ValueError("termination_shift must be a JSON object")
    dp1 = _config_value(config, ("dp1",), shifts.get("dp1", 0.0))
    dp2 = _config_value(config, ("dp2",), shifts.get("dp2", 0.0))
    return _finite_float(dp1, "dp1"), _finite_float(dp2, "dp2")


def _prepare_boundary(
    input_path: Path,
    config: Mapping[str, Any],
) -> tuple[Any, dict[str, Any], dict[str, Any]]:
    """Compute CSL, CNID, and shared settings once for a boundary."""

    np, Structure, get_primitive_hkl, rot, convert_vector_index, core = _backend()
    structure = Structure.from_file(input_path)

    axis = _config_value(config, ("rotation_axis", "axis"), None)
    if axis is None:
        raise ValueError("build mode requires rotation_axis, for example [0, 0, 1]")
    plane = _config_value(
        config,
        ("gb_plane", "grain_boundary_plane", "plane", "hkl"),
        None,
    )
    if plane is None:
        raise ValueError("build mode requires gb_plane, for example [1, 0, 0]")
    axis_input = axis
    plane_input_original = plane
    axis = three_indices(
        axis, "direction", config.get("rotation_axis_notation"), structure.lattice.matrix.T
    )
    plane = three_indices(
        plane, "plane", config.get("gb_plane_notation"), structure.lattice.matrix.T
    )
    target_sigma = _positive_int(
        _config_value(config, ("sigma", "target_sigma"), None),
        "sigma",
    )
    max_sigma = _positive_int(
        _config_value(config, ("max_sigma", "maxsigma"), target_sigma),
        "max_sigma",
    )
    angle_deg = float(
        _config_value(
            config,
            ("rotation_angle_deg", "angle_deg", "rotation_angle", "theta"),
            None,
        )
    )
    if not math.isfinite(angle_deg):
        raise ValueError("rotation_angle_deg must be finite")
    tolerances = _tolerances(config)
    csl_mode, exact = _csl_mode(config)

    interface, row = _search_angle(
        structure,
        axis,
        angle_deg,
        max_sigma,
        tolerances,
        exact,
        np,
        rot,
        convert_vector_index,
        core,
    )
    if row["status"] != "found":
        raise RuntimeError(
            f"no CSL candidate found for angle {angle_deg:g} degrees; "
            "increase max_sigma or the matching tolerances"
        )
    if row["sigma_1"] != target_sigma or row["sigma_2"] != target_sigma:
        raise RuntimeError(
            f"requested Sigma-{target_sigma}, but InterfaceMaster selected "
            f"Sigma ({row['sigma_1']}, {row['sigma_2']}); adjust tolerances "
            "or the requested angle"
        )

    plane_input = _finite_vector(plane, "gb_plane", 3, np)
    if np.linalg.norm(plane_input) == 0:
        raise ValueError("gb_plane must not be the zero vector")
    base = core(structure, structure, verbose=False)
    plane_primitive = get_primitive_hkl(
        plane_input,
        base.conv_lattice_1,
        base.lattice_1,
        tolerances["integer"],
    )

    normal_ortho = bool(config.get("normal_ortho", True))
    plane_ortho = bool(config.get("plane_ortho", False))
    superlattice_lim = _positive_int(
        config.get("superlattice_lim", config.get("lim", 50)),
        "superlattice_lim",
    )
    interface.compute_bicrystal(
        plane_primitive,
        lim=superlattice_lim,
        normal_ortho=normal_ortho,
        plane_ortho=plane_ortho,
        tol_ortho=float(config.get("tol_ortho", tolerances["integer"])),
        tol_integer=tolerances["integer"],
    )

    translation = _finite_vector(
        config.get("translation", config.get("dydz", [0.0, 0.0, 0.0])),
        "translation",
        3,
        np,
        allow_two=True,
    )
    normal_distance = _nonnegative_float(
        config.get("normal_distance", config.get("dx", 0.0)),
        "normal_distance",
    )
    dp1, dp2 = _termination_values(config)

    repeat_1 = _repeat_vector(
        config.get(
            "xyz_1",
            config.get("normal_repeat", config.get("repeat", [1, 1, 1])),
        ),
        "xyz_1",
    )
    repeat_2 = _repeat_vector(
        config.get(
            "xyz_2",
            config.get("normal_repeat", config.get("repeat", [1, 1, 1])),
        ),
        "xyz_2",
    )
    vacuum = _nonnegative_float(
        config.get("vacuum", config.get("vx", 0.0)),
        "vacuum",
    )
    kt_i = bool(config.get("KTI", config.get("kti", False)))
    mirror = bool(config.get("mirror", False))

    settings = {
        "dydz": translation,
        "dx": normal_distance,
        "dp1": dp1,
        "dp2": dp2,
        "xyz_1": repeat_1,
        "xyz_2": repeat_2,
        "vx": vacuum,
        "mirror": mirror,
        "KTI": kt_i,
    }
    report = {
        "input_structure": str(input_path),
        "rotation_axis": np.asarray(axis, dtype=float).tolist(),
        "rotation_axis_input": list(axis_input),
        "rotation_axis_notation": index_notation(axis_input, config.get("rotation_axis_notation")),
        "gb_plane_input": list(plane_input_original),
        "gb_plane_notation": index_notation(plane_input_original, config.get("gb_plane_notation")),
        "rotation_angle_deg": angle_deg,
        "csl_mode": csl_mode,
        "interface_master_exact": exact,
        "rotation_matrix_used": np.asarray(interface.R, dtype=float).tolist(),
        "sigma": target_sigma,
        "sigma_1": row["sigma_1"],
        "sigma_2": row["sigma_2"],
        "gb_plane": np.asarray(plane, dtype=float).tolist(),
        "gb_plane_primitive": np.asarray(plane_primitive, dtype=int).tolist(),
        "translation_cartesian_A": translation.tolist(),
        "normal_distance_A": normal_distance,
        "termination_shift": {"dp1": dp1, "dp2": dp2},
        "deformation_max_abs": row["deformation_max_abs"],
        "U1": np.asarray(interface.U1, dtype=int).tolist(),
        "U2": np.asarray(interface.U2, dtype=int).tolist(),
    }
    # mirror=True replaces grain 2 with grain 1's reflected slab. Its actual
    # surface translation lattice is then grain 1's, independent of R/D.
    if mirror:
        from interfacemaster.cellcalc import get_pri_vec_inplane

        interface.CNID = get_pri_vec_inplane(plane_primitive, interface.lattice_1)
    return interface, settings, report


def cnid_grid(config: Mapping[str, Any]) -> tuple[int, int]:
    """Read positive CNID subdivisions; each half-open grid has n1*n2 frames."""

    value = config.get("cnid_grid", [4, 4])
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ValueError("cnid_grid must contain [u1_subdivisions, u2_subdivisions]")
    return (
        _positive_int(value[0], "cnid_grid u1"),
        _positive_int(value[1], "cnid_grid u2"),
    )


def build_boundary_frames(
    input_path: Path,
    config: Mapping[str, Any],
    *,
    progress: Callable[[int, int], None] | None = None,
    cancelled: Callable[[], bool] | None = None,
) -> dict[str, Any]:
    """Sample a CNID grid using InterfaceMaster, returning independent frames.

    Frame i*n2+j has t=(i/n1)*u1+(j/n2)*u2 in the original Cartesian frame.
    The vectors are passed to get_bicrystal before its orientation adjustment.
    All frames use the same terminations and normal spacing. Legacy Cartesian
    translation keys remain supported by build_boundary but are ignored here.
    """

    subdivisions_1, subdivisions_2 = cnid_grid(config)
    input_path = input_path.expanduser().resolve()
    settings_config = dict(config)
    settings_config.pop("translation", None)
    settings_config.pop("dydz", None)
    np, _, _, _, _, _ = _backend()
    frames = []
    metadata = []
    total = subdivisions_1 * subdivisions_2
    if progress is not None:
        progress(0, total)
    with tempfile.TemporaryDirectory(prefix="cckit_cnid_") as work_dir:
        old_cwd = Path.cwd()
        try:
            os.chdir(work_dir)
            if cancelled is not None and cancelled():
                raise RuntimeError("CNID sampling was cancelled")
            interface, settings, report = _prepare_boundary(input_path, settings_config)
            vectors = np.asarray(interface.CNID, dtype=float)
            if vectors.shape != (3, 2) or not np.all(np.isfinite(vectors)):
                raise RuntimeError("InterfaceMaster did not return two finite CNID vectors")
            if np.linalg.norm(np.cross(vectors[:, 0], vectors[:, 1])) < 1e-12:
                raise RuntimeError("InterfaceMaster returned a degenerate CNID cell")
            for i in range(subdivisions_1):
                for j in range(subdivisions_2):
                    if cancelled is not None and cancelled():
                        raise RuntimeError("CNID sampling was cancelled")
                    alpha = i / subdivisions_1
                    beta = j / subdivisions_2
                    translation = vectors @ np.asarray([alpha, beta])
                    settings["dydz"] = translation
                    # One scratch POSCAR is overwritten; the GUI receives
                    # copies, not references to InterfaceMaster's mutable state.
                    interface.get_bicrystal(
                        **settings,
                        filename=str(Path(work_dir) / "POSCAR_frame"),
                        filetype="VASP",
                    )
                    frame = interface.bicrystal_structure.copy()
                    frames.append(frame)
                    metadata.append({
                        "frame": len(frames) - 1,
                        "cnid_i": i,
                        "cnid_j": j,
                        "cnid_n1": subdivisions_1,
                        "cnid_n2": subdivisions_2,
                        "rbt_alpha": alpha,
                        "rbt_beta": beta,
                        "rbt_translation_A": translation.tolist(),
                        "cnid_u1_A": vectors[:, 0].tolist(),
                        "cnid_u2_A": vectors[:, 1].tolist(),
                        "termination_dp1_A": settings["dp1"],
                        "termination_dp2_A": settings["dp2"],
                    })
                    if progress is not None:
                        progress(len(frames), total)
        finally:
            os.chdir(old_cwd)
    return {
        "frames": frames,
        "metadata": metadata,
        "cnid_vectors_cartesian_A": vectors.tolist(),
        "report": report,
    }


def build_boundary(
    input_path: Path,
    config: Mapping[str, Any],
    output_path: Path,
    report_path: Path | None,
) -> Path:
    """Build one boundary; preserve the public CLI's single-POSCAR semantics."""

    interface, settings, report = _prepare_boundary(input_path, config)
    _, Structure, _, _, _, _ = _backend()
    # InterfaceMaster writes four grain CIFs beside its temporary POSCAR.
    # Keep both the pre-termination cells and the slabs after dp1/dp2 so the
    # two cut surfaces can be inspected independently of the combined model.
    grain_cifs = {
        name: output_path.with_name(f"{output_path.name}_{name}")
        for name in (
            "min_sub.cif",
            "min_film.cif",
            "super_sub.cif",
            "super_film.cif",
        )
    }
    with tempfile.TemporaryDirectory(prefix="cckit_gb_") as work_dir:
        old_cwd = Path.cwd()
        try:
            os.chdir(work_dir)
            interface.get_bicrystal(
                **settings,
                filename=str(output_path),
                filetype="VASP",
            )
        finally:
            os.chdir(old_cwd)

        missing = [name for name in grain_cifs if not (Path(work_dir) / name).is_file()]
        if missing:
            raise RuntimeError(
                "InterfaceMaster did not create the grain CIF files: "
                + ", ".join(missing)
            )
        for name, destination in grain_cifs.items():
            shutil.copy2(Path(work_dir) / name, destination)

    if not output_path.is_file():
        raise RuntimeError(f"InterfaceMaster did not create {output_path}")

    final_structure = Structure.from_file(output_path)
    report.update({
        "output_poscar": str(output_path),
        "grain_cifs": {name: str(path) for name, path in grain_cifs.items()},
        "atom_count": len(final_structure),
    })
    if report_path is not None:
        report_path.write_text(
            json.dumps(report, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )

    print(
        f"Built Sigma-{report['sigma']} "
        f"{report['rotation_angle_deg']:g} degree symmetric tilt "
        f"boundary ({len(final_structure)} atoms)."
    )
    print(f"POSCAR: {output_path}")
    for name, path in grain_cifs.items():
        print(f"{name}: {path}")
    if report_path is not None:
        print(f"Report: {report_path}")
    return output_path


def _filename_number(value: float) -> str:
    """Format a numeric filename component without a decimal point."""

    if value == 0:
        return "0"
    return f"{value:g}".replace(".", "_")


def _default_output_name(mode: str, config: Mapping[str, Any]) -> str:
    """Return a deterministic default output filename."""

    if mode == "scan":
        return str(config.get("scan_output", "GB_scan_results.csv"))
    sigma = config.get("sigma", config.get("target_sigma", "unknown"))
    angle = config.get(
        "rotation_angle_deg",
        config.get(
            "angle_deg",
            config.get("rotation_angle", config.get("theta", "unknown")),
        ),
    )
    try:
        angle_value = float(angle)
    except (TypeError, ValueError):
        angle_name = str(angle)
    else:
        angle_name = (
            str(int(round(angle_value)))
            if math.isfinite(angle_value)
            else "unknown"
        )
    dp1, dp2 = _termination_values(config)
    return (
        f"POSCAR_Sigma{sigma}_{angle_name}deg_"
        f"dp1_{_filename_number(dp1)}_dp2_{_filename_number(dp2)}_GB"
    )


def _parse_args() -> argparse.Namespace:
    """Parse the public command-line interface."""

    parser = argparse.ArgumentParser(
        description=(
            "Scan CSL candidates or build a symmetric-tilt grain boundary "
            "with InterfaceMaster and gbin.json."
        )
    )
    parser.add_argument("input", help="input crystal structure file")
    parser.add_argument(
        "--config",
        default=None,
        help=(
            "configuration JSON (default: gbin.json beside the input or in "
            "the current directory)"
        ),
    )
    parser.add_argument(
        "--mode",
        choices=("scan", "build"),
        default=None,
        help="override gbin.json mode",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help=(
            "output CSV or POSCAR; a bare filename is written beside the "
            "input (build names receive a POSCAR_ prefix)"
        ),
    )
    parser.add_argument(
        "--report",
        default=None,
        help="optional JSON report for build mode",
    )
    return parser.parse_args()


def main() -> None:
    """Run scan or build mode and convert expected failures to CLI errors."""

    args = _parse_args()
    try:
        input_path = _resolve_input(args.input)
        config_path = _resolve_config(args.config, input_path)
        config = _read_json(config_path)
        mode = args.mode or config.get("mode", "scan")
        if mode not in {"scan", "build"}:
            raise ValueError("mode must be either 'scan' or 'build'")

        configured_output = args.output or _config_value(
            config,
            (
                "output",
                "output_file",
                "scan_output" if mode == "scan" else "poscar_output",
            ),
            None,
        )
        output = _resolve_output(
            configured_output,
            input_path,
            _default_output_name(mode, config),
        )
        if mode == "build":
            output = _ensure_poscar_prefix(output)
            report_value = args.report or config.get("report")
            report = (
                _resolve_output(report_value, input_path, f"{output.name}.json")
                if report_value
                else None
            )
            build_boundary(input_path, config, output, report)
        else:
            scan_candidates(input_path, config, output)
    except (FileNotFoundError, RuntimeError, ValueError, TypeError) as exc:
        raise SystemExit(f"model-gb: error: {exc}") from exc


if __name__ == "__main__":
    main()

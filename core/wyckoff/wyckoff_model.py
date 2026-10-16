# CCkit 功能说明：根据空间群、晶格参数和 Wyckoff 代表坐标生成晶体。
"""Deterministic Wyckoff expansion in conventional spglib settings.

Lengths are in Angstrom, angles in degrees, and site coordinates fractional.
Monoclinic cells use unique axis b; trigonal R groups use hexagonal axes.
"""

from __future__ import annotations

import math
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable

import numpy as np
from pymatgen.core import Element, Lattice, Structure


def coordinate_value(value: Any) -> float:
    """Read a decimal or a simple fraction without evaluating Python code."""

    try:
        number = float(Fraction(str(value).strip()))
    except (ValueError, ZeroDivisionError, TypeError, OverflowError) as exc:
        raise ValueError(f"Invalid number {value!r}; use a decimal or fraction such as 1/3") from exc
    if not math.isfinite(number):
        raise ValueError("Coordinates and lattice parameters must be finite")
    return number


@lru_cache(maxsize=230)
def space_group(number: int) -> Any:
    """Return space-group data with one consistent conventional setting."""

    if isinstance(number, bool) or not isinstance(number, int) or not 1 <= number <= 230:
        raise ValueError("Space group number must be an integer from 1 to 230")
    try:
        from pyxtal.symmetry import Group
    except ImportError as exc:
        raise RuntimeError(
            "Wyckoff model requires PyXtal. Install requirements-gb-viewer.txt "
            "in the Python environment running this viewer."
        ) from exc
    return Group(number, style="spglib")


def lattice_fields(number: int) -> tuple[str, ...]:
    """List the independent conventional-cell parameters for a space group."""

    space_group(number)
    if number <= 2:
        return ("a", "b", "c", "alpha", "beta", "gamma")
    if number <= 15:
        return ("a", "b", "c", "beta")
    if number <= 74:
        return ("a", "b", "c")
    if number <= 194:
        return ("a", "c")
    return ("a",)


def lattice_from_config(number: int, values: dict[str, Any]) -> Lattice:
    """Enforce crystal-system constraints and reject degenerate cells."""

    parameters = {key: coordinate_value(values[key]) for key in lattice_fields(number)}
    for key, value in parameters.items():
        if key in {"a", "b", "c"} and value <= 0:
            raise ValueError(f"Lattice constant {key} must be positive")
        if key in {"alpha", "beta", "gamma"} and not 0 < value < 180:
            raise ValueError(f"Angle {key} must be between 0 and 180 degrees")
    a = parameters["a"]
    b = parameters.get("b", a)
    c = parameters.get("c", a)
    alpha = parameters.get("alpha", 90.0)
    beta = parameters.get("beta", 90.0)
    gamma = parameters.get("gamma", 120.0 if 143 <= number <= 194 else 90.0)
    cosines = np.cos(np.radians([alpha, beta, gamma]))
    metric_factor = 1 + 2 * np.prod(cosines) - np.sum(cosines ** 2)
    if metric_factor <= 1e-10:
        raise ValueError("Lattice angles do not define a cell with positive volume")
    return Lattice.from_parameters(a, b, c, alpha, beta, gamma)


def default_config() -> dict[str, Any]:
    """Provide the hcp Zr example used by the experimental Wyckoff script."""

    return {
        "mode": "build", "space_group": 194,
        "lattice": {"a": 3.231, "c": 5.147},
        "sites": [{
            "element": "Zr", "wyckoff": "2c", "coords": ["1/3", "2/3", "1/4"],
            "coordinate_mode": "equivalent",
        }],
    }


def parameters(config: dict[str, Any]) -> dict[str, Any]:
    """Validate and normalize the card configuration before building."""

    number = config.get("space_group")
    group = space_group(number)
    lattice_values = config.get("lattice")
    if not isinstance(lattice_values, dict):
        raise ValueError("Lattice parameters must be an object")
    missing = set(lattice_fields(number)) - lattice_values.keys()
    if missing:
        raise ValueError(f"Missing lattice parameters: {', '.join(sorted(missing))}")
    lattice_from_config(number, lattice_values)
    sites = config.get("sites")
    if not isinstance(sites, list) or not sites:
        raise ValueError("Add at least one Wyckoff position")
    labels = {position.get_label() for position in group}
    normalized = []
    for index, site in enumerate(sites, 1):
        if not isinstance(site, dict):
            raise ValueError(f"Site {index}: expected element, Wyckoff position and coordinates")
        element = str(site.get("element", "")).strip()
        if not Element.is_valid_symbol(element):
            raise ValueError(f"Site {index}: invalid element {element!r} (examples: Zr, O, Si)")
        label = site.get("wyckoff")
        if label not in labels:
            raise ValueError(f"Site {index}: {label!r} is not a Wyckoff position in group {number}")
        coords = site.get("coords")
        if not isinstance(coords, (list, tuple)) or len(coords) != 3:
            raise ValueError(f"Site {index}: provide fractional x, y and z")
        coordinate_mode = site.get("coordinate_mode", "manual")
        if coordinate_mode not in {"equivalent", "manual"}:
            raise ValueError(f"Site {index}: coordinate_mode must be 'equivalent' or 'manual'")
        normalized.append({
            "element": element, "wyckoff": label,
            "coords": [coordinate_value(value) for value in coords],
            "coordinate_mode": coordinate_mode,
        })
    return {
        "mode": "build", "space_group": number,
        "lattice": {key: coordinate_value(lattice_values[key]) for key in lattice_fields(number)},
        "sites": normalized,
    }


def _unique_coordinates(coords: np.ndarray) -> np.ndarray:
    """Deduplicate fractional coordinates using periodic boundaries."""

    unique = []
    for point in np.asarray(coords) % 1.0:
        if unique:
            delta = np.asarray(unique) - point
            delta -= np.rint(delta)
            if np.any(np.max(np.abs(delta), axis=1) < 1e-5):
                continue
        unique.append(point)
    return np.asarray(unique)


def position_options(
    number: int, label: str, coords: list[Any] | None = None,
    *, free_values: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """List a complete Wyckoff orbit and its symbolic coordinate expressions.

    Free parameters describe a continuum of different orbits. With no supplied
    coordinates, a deterministic example is used, rather than an enumeration
    of every possible free-parameter value. Any equivalent orbit member may
    be supplied as the representative coordinate. ``free_values`` instead
    supplies the x/y/z variables used in the symbolic Wyckoff operators.
    """

    group = space_group(number)
    position = next((position for position in group if position.get_label() == label), None)
    if position is None:
        raise ValueError(f"{label!r} is not a Wyckoff position in group {number}")
    operation = position.ops[0]
    axes = [axis for axis in range(3) if any(
        np.any(np.abs(item.rotation_matrix[:, axis]) > 1e-12) for item in position.ops
    )]
    names = ["xyz"[axis] for axis in axes]
    raw_parameters = np.asarray([0.137, 0.271, 0.389])
    if free_values is not None:
        if coords is not None:
            raise ValueError("Supply either representative coordinates or free parameters")
        if set(free_values) != set(names):
            raise ValueError(f"Provide exactly these Wyckoff parameters: {', '.join(names) or 'none'}")
        raw_parameters = np.zeros(3)
        for axis, name in zip(axes, names):
            raw_parameters[axis] = coordinate_value(free_values[name])
    if coords is None:
        point = operation.operate(raw_parameters) % 1.0
    else:
        if len(coords) != 3:
            raise ValueError("Provide fractional x, y and z")
        point = np.asarray([coordinate_value(value) for value in coords]) % 1.0
    representative = position.search_generator(
        point, ops=group[0].ops, tol=1e-5, symmetrize=True,
    )
    if representative is None:
        raise ValueError(
            f"Coordinates do not match {label}. "
            f"Representative constraint: {position.ops[0].as_xyz_str()}"
        )
    if coords is not None:
        # Recover values of the variables in the symbolic operators, rather
        # than confusing them with the components of an equivalent point.
        raw_parameters = np.zeros(3)
        if axes:
            values, *_ = np.linalg.lstsq(
                operation.rotation_matrix[:, axes],
                representative - operation.translation_vector, rcond=None,
            )
            raw_parameters[axes] = values
    orbit = _unique_coordinates(position.apply_ops(raw_parameters))
    full_orbit = _unique_coordinates(group[0].apply_ops(representative))
    if len(orbit) != position.multiplicity or len(full_orbit) != position.multiplicity:
        raise ValueError(
            f"Coordinates have a different multiplicity from {label}; "
            "choose the matching special Wyckoff position"
        )
    return {
        "label": label,
        "multiplicity": position.multiplicity,
        "free_parameters": int(position.get_dof()),
        "parameter_names": names,
        "parameter_values": {name: float(raw_parameters[axis]) for axis, name in zip(axes, names)},
        "constraint": position.ops[0].as_xyz_str(),
        "equivalents": [{
            "expression": operation.as_xyz_str(), "coords": coordinate.tolist(),
        } for operation, coordinate in zip(position.ops, orbit)],
    }


def build_wyckoff_frames(
    _input_path: Path, config: dict[str, Any], *,
    progress: Callable[[int, int], None] | None = None,
    cancelled: Callable[[], bool] | None = None,
) -> dict[str, Any]:
    """Expand each chosen orbit into one fully occupied conventional cell."""

    values = parameters(config)
    group = space_group(values["space_group"])
    lattice = lattice_from_config(group.number, values["lattice"])
    species: list[str] = []
    coordinates: list[Any] = []
    for index, site in enumerate(values["sites"], 1):
        if cancelled is not None and cancelled():
            raise RuntimeError("Wyckoff generation cancelled")
        try:
            options = position_options(group.number, site["wyckoff"], site["coords"])
        except ValueError as exc:
            raise ValueError(f"Site {index}: {exc}") from exc
        for equivalent in options["equivalents"]:
            point = np.asarray(equivalent["coords"])
            if coordinates:
                delta = np.asarray(coordinates) - point
                delta -= np.rint(delta)
                if np.any(np.max(np.abs(delta), axis=1) < 1e-5):
                    raise ValueError(f"Site {index}: overlaps another Wyckoff orbit; remove duplicate sites")
            coordinates.append(point)
            species.append(site["element"])
        if progress is not None:
            progress(index, len(values["sites"]))
    structure = Structure(lattice, species, coordinates)
    return {
        "frames": [structure],
        "metadata": [{"frame": 0, "space_group": group.number, "interface_plane_fractions": []}],
        "report": {"space_group": group.number, "symbol": group.symbol, "atom_count": len(structure)},
    }

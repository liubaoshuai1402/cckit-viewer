"""Cut bulk-periodic slabs with explicitly prescribed in-plane crystal axes."""

from __future__ import annotations

import copy
import math
from pathlib import Path
from typing import Any, Callable, Mapping

from core.bicrystal.bicrystal_backend import index_notation, index_values, three_indices


def default_config() -> dict[str, Any]:
    """Return an independent (001) slab with input-cell [100]/[010] axes."""

    return {
        "mode": "build", "bulk_structure": "",
        "plane": [0, 0, 1], "plane_notation": "three",
        "direction_a": [1, 0, 0], "direction_a_notation": "three",
        "direction_b": [0, 1, 0], "direction_b_notation": "three",
        "repeat_a": 1, "repeat_b": 1,
        "thickness": 15.0, "vacuum": 15.0,
        "termination_mode": "layer", "termination_index": 0, "cut_phase": 0.0,
        "layer_tolerance": 0.1, "max_atoms": 50000,
    }


def parameters(config: Mapping[str, Any]) -> dict[str, Any]:
    """Validate JSON controls without loading a crystal or allocating a slab."""

    if not isinstance(config, Mapping):
        raise ValueError("surface.json must contain a JSON object")
    values = default_config()
    values.update(copy.deepcopy(dict(config)))
    if values["mode"] != "build":
        raise ValueError("Surface Builder supports build mode only")
    if not isinstance(values["bulk_structure"], str):
        raise ValueError("bulk_structure must be a file path")
    for key in ("plane", "direction_a", "direction_b"):
        notation = index_notation(values[key], config.get(f"{key}_notation"))
        indices = index_values(values[key], "plane" if key == "plane" else "direction", notation)
        if any(item != round(item) or abs(item) > 2000000 for item in indices):
            raise ValueError(f"{key} requires integer indices between -2000000 and 2000000")
        values[key], values[f"{key}_notation"] = [int(item) for item in indices], notation
    for key, lower, upper in (
        ("repeat_a", 1, 100), ("repeat_b", 1, 100),
        ("termination_index", 0, 9999), ("max_atoms", 1, 200000),
    ):
        value = values[key]
        if isinstance(value, bool) or not isinstance(value, int) or not lower <= value <= upper:
            raise ValueError(f"{key} must be an integer between {lower} and {upper}")
    for key, lower, upper in (
        ("thickness", 0.001, 10000), ("vacuum", 0.001, 10000),
        ("cut_phase", 0, 1), ("layer_tolerance", 0.000001, 5),
    ):
        raw = values[key]
        if isinstance(raw, bool):
            raise ValueError(f"{key} must be a finite number")
        try:
            value = float(raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{key} must be a finite number") from exc
        if not math.isfinite(value) or not lower <= value <= upper:
            raise ValueError(f"{key} must be between {lower} and {upper}")
        values[key] = value
    values["cut_phase"] %= 1.0
    if values["termination_mode"] not in {"layer", "custom"}:
        raise ValueError("termination_mode must be layer or custom")
    return values


def _check_cancelled(cancelled: Callable[[], bool] | None) -> None:
    if cancelled is not None and cancelled():
        raise RuntimeError("Surface construction cancelled")


def _primitive_indices(values: list[float]) -> list[int]:
    indices = [int(round(item)) for item in values]
    divisor = math.gcd(*indices)
    return [item // divisor for item in indices]


def _bezout(first: int, second: int) -> tuple[int, int, int]:
    """Return g, x, y with first*x + second*y = gcd(first, second)."""

    old_r, remainder = abs(first), abs(second)
    old_x, x, old_y, y = 1, 0, 0, 1
    while remainder:
        quotient = old_r // remainder
        old_r, remainder = remainder, old_r - quotient * remainder
        old_x, x = x, old_x - quotient * x
        old_y, y = y, old_y - quotient * y
    return old_r, old_x * (1 if first >= 0 else -1), old_y * (1 if second >= 0 else -1)


def _oriented_cell(bulk: Any, values: Mapping[str, Any]) -> tuple[Any, Any, Any, float]:
    """Construct an integer bulk supercell while retaining both requested axes."""

    import numpy as np

    converted = {}
    for key in ("plane", "direction_a", "direction_b"):
        converted[key] = _primitive_indices(three_indices(
            values[key], "plane" if key == "plane" else "direction",
            values[f"{key}_notation"], bulk.lattice.matrix.T,
        ))
    h = converted["plane"]
    u, v = converted["direction_a"], converted["direction_b"]
    for name, vector in (("a", u), ("b", v)):
        if sum(first * second for first, second in zip(h, vector)) != 0:
            raise ValueError(f"Direction {name} is not in the selected plane: hU + kV + lW must be zero")
    cross = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
    if not any(cross):
        raise ValueError("The two in-plane directions must not be parallel")
    gcd_hk, x, y = _bezout(h[0], h[1])
    _gcd_all, scale, z = _bezout(gcd_hk, h[2])
    w = np.array([x * scale, y * scale, z], dtype=np.int64)
    u = np.array(u, dtype=np.int64) * values["repeat_a"]
    v = np.array(v, dtype=np.int64) * values["repeat_b"]
    matrix = np.vstack((u, v, w))
    # Integer volume is computed before pymatgen allocates supercell sites.
    multiplicity = sum(int(a) * int(b) for a, b in zip(np.cross(u, v), w))
    if multiplicity == 0:
        raise ValueError("The requested surface directions do not span a valid periodic cell")
    if abs(multiplicity) * len(bulk) > values["max_atoms"]:
        raise ValueError("The oriented cell exceeds Maximum atoms; use smaller direction indices or repeats")
    lattice = np.asarray(bulk.lattice.matrix)
    normal = np.linalg.solve(lattice, h)
    period = 1.0 / float(np.linalg.norm(normal))
    normal /= np.linalg.norm(normal)
    a, b = u @ lattice, v @ lattice
    if np.dot(np.cross(a, b), normal) <= 0:
        raise ValueError("The surface directions must form a right-handed pair with the plane normal; swap a/b or reverse one direction")
    # Shorten only the auxiliary stacking vector; the requested axes stay fixed.
    correction = np.rint(np.linalg.lstsq(np.column_stack((a, b)), w @ lattice, rcond=None)[0]).astype(np.int64)
    matrix[2] -= correction[0] * u + correction[1] * v
    # pymatgen enumerates a bounding box before retaining supercell points.
    # A highly skew integer matrix can make that box huge despite a small det.
    enumeration_box = [sum(abs(int(item)) for item in matrix[:, axis]) + 1 for axis in range(3)]
    if math.prod(enumeration_box) > 2000000:
        raise ValueError("These direction indices create an excessively skew supercell; choose shorter in-plane directions")
    oriented = bulk.copy()
    oriented.make_supercell(matrix)
    return oriented, matrix, normal, period


def _cut_choices(oriented: Any, normal: Any, period: float, tolerance: float) -> list[dict[str, Any]]:
    """List gaps between periodic atomic layers; no symmetry-equivalent cuts are dropped."""

    import numpy as np

    phases = (oriented.cart_coords @ normal / period) % 1.0
    phases[np.isclose(phases, 1, rtol=0, atol=1e-8)] = 0
    heights = np.unique(np.round(phases, 10))
    choices = []
    for index, height in enumerate(heights):
        next_height = heights[(index + 1) % len(heights)] + (1 if index == len(heights) - 1 else 0)
        if (next_height - height) * period <= tolerance:
            continue
        cut = float((height + next_height) / 2 % 1)
        bottom = float(next_height % 1)
        distance = (phases - bottom) % 1
        distance[np.isclose(distance, 1, rtol=0, atol=1e-8)] = 0
        species = sorted({str(site.specie) for site, delta in zip(oriented, distance)
                          if delta * period <= tolerance + 1e-8})
        choices.append({"phase": cut, "bottom_species": species})
    if len(choices) > 10000:
        raise ValueError("More than 10000 distinct cuts; reduce the input cell or increase Layer tolerance")
    return sorted(choices, key=lambda choice: choice["phase"])


def build_surface_frames(
    input_path: Path, config: Mapping[str, Any], *,
    progress: Callable[[int, int], None] | None = None,
    cancelled: Callable[[], bool] | None = None,
) -> dict[str, Any]:
    """Slice an unstrained periodic bulk, then add an exact atom-to-atom vacuum gap.

    Thickness is the normal width of a half-open cutting window [cut, cut+t).
    Its atomic extent is discrete and can be smaller. Vacuum is the normal
    distance between outermost atomic planes of neighboring periodic slabs.
    """

    import numpy as np
    from pymatgen.core import Structure

    values = parameters(config)
    _check_cancelled(cancelled)
    if progress:
        progress(0, 1)
    bulk = Structure.from_file(Path(input_path).expanduser().resolve())
    if not bulk.is_ordered or not len(bulk):
        raise ValueError("Choose a nonempty bulk structure with fully occupied atomic sites")
    oriented, matrix, normal, period = _oriented_cell(bulk, values)
    _check_cancelled(cancelled)
    choices = _cut_choices(oriented, normal, period, values["layer_tolerance"])
    if values["termination_mode"] == "layer":
        if not choices:
            raise ValueError("No interlayer gap exceeds Layer tolerance; reduce that tolerance or use a custom cut")
        if values["termination_index"] >= len(choices):
            raise ValueError(f"Only {len(choices)} cuts are available; choose an index from 1 to {len(choices)}")
        phase = choices[values["termination_index"]]["phase"]
    else:
        phase = values["cut_phase"]
    cycles = math.ceil(values["thickness"] / period) + 1
    if cycles * len(oriented) > 1000000:
        raise ValueError("The requested slab exceeds the construction limit; reduce thickness or in-plane repeats")
    start, stop = phase * period, phase * period + values["thickness"]
    coordinates, species, source_indices = [], [], []
    # Wrapping fixes numerical points at a repeat boundary, retaining complete layers.
    fractional = oriented.frac_coords % 1.0
    fractional[np.isclose(fractional, 1, rtol=0, atol=1e-8)] = 0
    base = fractional @ oriented.lattice.matrix
    for cycle in range(cycles + 1):
        _check_cancelled(cancelled)
        positions = base + cycle * oriented.lattice.matrix[2]
        heights = positions @ normal
        selected = np.flatnonzero((heights >= start - 1e-8) & (heights < stop - 1e-8))
        if len(species) + len(selected) > values["max_atoms"]:
            raise ValueError("The surface exceeds Maximum atoms; reduce thickness or in-plane repeats")
        coordinates.extend(positions[selected])
        species.extend(oriented[int(index)].specie for index in selected)
        source_indices.extend(int(index) for index in selected)
    if not species:
        raise ValueError("The cutting window contains no atoms; increase thickness or change the cut")
    coords = np.asarray(coordinates)
    heights = coords @ normal
    bottom, top = float(heights.min()), float(heights.max())
    extent = top - bottom
    height = extent + values["vacuum"]
    # Rigid rotation: a along x, b in xy, surface normal along +z.
    a, b = oriented.lattice.matrix[:2]
    x_axis = a / np.linalg.norm(a)
    frame = np.vstack((x_axis, np.cross(normal, x_axis), normal))
    coords = coords @ frame.T
    coords[:, 2] += values["vacuum"] / 2 - bottom
    lattice = np.vstack((a @ frame.T, b @ frame.T, [0, 0, height]))
    frac = coords @ np.linalg.inv(lattice)
    frac[:, :2] %= 1.0
    properties = {key: [items[index] for index in source_indices]
                  for key, items in oriented.site_properties.items()}
    structure = Structure(lattice, species, frac, site_properties=properties)
    top_species = sorted({str(atom) for atom, value in zip(species, heights)
                          if top - value <= values["layer_tolerance"] + 1e-8})
    bottom_species = sorted({str(atom) for atom, value in zip(species, heights)
                             if value - bottom <= values["layer_tolerance"] + 1e-8})
    planes = [values["vacuum"] / 2 / height, (values["vacuum"] / 2 + extent) / height]
    _check_cancelled(cancelled)
    if progress:
        progress(1, 1)
    report = {
        "atom_count": len(structure), "area_A2": float(np.linalg.norm(np.cross(a, b))),
        "cut_window_A": values["thickness"], "atomic_thickness_A": extent,
        "vacuum_A": values["vacuum"], "normal_period_A": period,
        "cut_phase": phase, "cut_choices": choices,
        "top_species": top_species, "bottom_species": bottom_species,
        "oriented_matrix": matrix.tolist(), "in_plane_vectors_A": lattice[:2].tolist(),
        "config": copy.deepcopy(values),
    }
    return {
        "frames": [structure],
        "metadata": [{"frame": 0, "feature": "surface", "interface_plane_fractions": [],
                      "surface_plane_fractions": planes, "cut_phase": phase}],
        "report": report,
    }


def _default_output_name(mode: str, config: Mapping[str, Any]) -> str:
    """Name exported slabs without referring to another modeling feature."""

    return "POSCAR_Surface"

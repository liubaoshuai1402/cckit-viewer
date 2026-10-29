"""Zur--McGill surface matching and coherent atomic interface construction.

A small cache reuses expensive searches when only slab geometry or the selected
match changes. Cached builders are never modified during atomic builds.
"""

from __future__ import annotations

import copy
import json
import math
from collections import OrderedDict
from pathlib import Path
from threading import Lock
from typing import Any, Callable, Iterator, Mapping

from core.interface.interface_backend import SEARCH_KEYS, frame_identity, shift_grid


_SEARCH_CACHE: OrderedDict[str, tuple[Any, dict[str, Any]]] = OrderedDict()
_CACHE_LOCK = Lock()


def _check_cancelled(cancelled: Callable[[], bool] | None) -> None:
    if cancelled is not None and cancelled():
        raise RuntimeError("Interface operation cancelled")


def _paths(input_path: Path, config: Mapping[str, Any]) -> tuple[Path, Path]:
    substrate = input_path.expanduser().resolve()
    if not config["film_structure"].strip():
        raise ValueError("Choose a film structure before searching")
    film = Path(config["film_structure"]).expanduser()
    if not film.is_absolute():
        film = substrate.parent / film
    film = film.resolve()
    for label, path in (("substrate", substrate), ("film", film)):
        if not path.is_file():
            raise FileNotFoundError(f"{label} structure was not found: {path}")
    return substrate, film


def _signature(input_path: Path, config: Mapping[str, Any]) -> str:
    paths = _paths(input_path, config)
    files = [(str(p), p.stat().st_mtime_ns, p.stat().st_size) for p in paths]
    return json.dumps([files, {key: config[key] for key in SEARCH_KEYS}], sort_keys=True)


def _miller(structure: Any, config: Mapping[str, Any], side: str) -> tuple[int, int, int]:
    indices = config[f"{side}_plane"]
    if len(indices) == 4:
        lattice = structure.lattice
        if not (
            math.isclose(lattice.a, lattice.b, rel_tol=1e-5)
            and abs(lattice.alpha - 90) < 0.01 and abs(lattice.beta - 90) < 0.01
            and abs(lattice.gamma - 120) < 0.01
        ):
            raise ValueError(f"{side}: four-index planes require a standard hexagonal input cell")
        indices = [indices[0], indices[1], indices[3]]
    return tuple(indices)


def _match_record(match: Any, index: int) -> dict[str, Any]:
    import numpy as np
    from pymatgen.analysis.elasticity.strain import Deformation
    from pymatgen.analysis.interfaces.coherent_interfaces import get_2d_transform

    film = np.asarray(match.film_sl_vectors)
    substrate = np.asarray(match.substrate_sl_vectors)
    strain = Deformation(match.match_transformation).green_lagrange_strain

    def angle(vectors: Any) -> float:
        return float(np.arctan2(np.linalg.norm(np.cross(*vectors)), np.dot(*vectors)))

    film_angle, substrate_angle = angle(film), angle(substrate)
    return {
        "index": index,
        "film_area_A2": float(np.linalg.norm(np.cross(*film))),
        "substrate_area_A2": float(np.linalg.norm(np.cross(*substrate))),
        "length_mismatch": (np.linalg.norm(substrate, axis=1) / np.linalg.norm(film, axis=1) - 1).tolist(),
        "angle_mismatch": substrate_angle / film_angle - 1,
        "angle_difference_deg": math.degrees(substrate_angle - film_angle),
        "von_mises_strain": float(strain.von_mises_strain),
        "strain": np.asarray(strain).tolist(),
        "deformation": np.asarray(match.match_transformation).tolist(),
        "film_vectors_A": film.tolist(), "substrate_vectors_A": substrate.tolist(),
        # ZSL's HNF matrices precede basis reduction. Include reduction here,
        # matching the integer transforms used by CoherentInterfaceBuilder.
        "film_supercell_matrix": np.rint(get_2d_transform(match.film_vectors, film)).astype(int).tolist(),
        "substrate_supercell_matrix": np.rint(get_2d_transform(match.substrate_vectors, substrate)).astype(int).tolist(),
    }


def _prepare(
    input_path: Path, config: Mapping[str, Any], cancelled: Callable[[], bool] | None,
) -> tuple[Any, dict[str, Any]]:
    from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder
    from pymatgen.analysis.interfaces.zsl import ZSLGenerator
    from pymatgen.core import Structure

    _check_cancelled(cancelled)
    signature = _signature(input_path, config)
    with _CACHE_LOCK:
        cached = _SEARCH_CACHE.get(signature)
        if cached is not None:
            _SEARCH_CACHE.move_to_end(signature)
            return cached

    class LimitedZSL(ZSLGenerator):
        """Bound retained results and check cancellation between search sets."""

        truncated = False

        def generate_sl_transformation_sets(self, film_area: Any, substrate_area: Any) -> Iterator[Any]:
            for transforms in super().generate_sl_transformation_sets(film_area, substrate_area):
                _check_cancelled(cancelled)
                yield transforms

        def __call__(self, film_vectors: Any, substrate_vectors: Any, lowest: bool = False) -> Iterator[Any]:
            matches = super().__call__(film_vectors, substrate_vectors, lowest=config["lowest"])
            for index, match in enumerate(matches):
                _check_cancelled(cancelled)
                if index >= config["max_matches"]:
                    self.truncated = True
                    break
                yield match

    substrate_path, film_path = _paths(input_path, config)
    substrate, film = Structure.from_file(substrate_path), Structure.from_file(film_path)
    substrate_plane = _miller(substrate, config, "substrate")
    film_plane = _miller(film, config, "film")
    generator = LimitedZSL(**{key: config[key] for key in (
        "max_area_ratio_tol", "max_area", "max_length_tol", "max_angle_tol", "bidirectional",
    )})
    builder = CoherentInterfaceBuilder(
        substrate_structure=substrate, film_structure=film,
        film_miller=film_plane, substrate_miller=substrate_plane, zslgen=generator,
        termination_ftol=(config["film_termination_ftol"], config["substrate_termination_ftol"]),
        label_index=config["label_index"], filter_out_sym_slabs=config["filter_out_sym_slabs"],
    )
    _check_cancelled(cancelled)
    if not builder.zsl_matches:
        raise ValueError("No lattice matches found. Increase the search area or relax matching tolerances.")
    if not builder.terminations:
        raise ValueError("No surface terminations found for these planes and tolerances")
    records = [_match_record(match, i) for i, match in enumerate(builder.zsl_matches)]
    ordered = sorted(zip(builder.zsl_matches, records), key=lambda pair: (
        pair[1]["film_area_A2"], pair[1]["von_mises_strain"],
    ))
    builder.zsl_matches = [match for match, _record in ordered]
    records = []
    for i, (_match, record) in enumerate(ordered):
        record["index"] = i
        records.append(record)
    result = {
        "matches": records,
        "terminations": [list(pair) for pair in builder.terminations],
        "truncated": generator.truncated,
        "substrate_miller": list(substrate_plane), "film_miller": list(film_plane),
    }
    _check_cancelled(cancelled)
    # Cached native generators must not retain a finished Qt worker callback.
    cancelled = None
    with _CACHE_LOCK:
        _SEARCH_CACHE[signature] = (builder, result)
        _SEARCH_CACHE.move_to_end(signature)
        while len(_SEARCH_CACHE) > 2:
            _SEARCH_CACHE.popitem(last=False)
    return builder, result


def search_matches(
    input_path: Path, config: Mapping[str, Any], *, cancelled: Callable[[], bool] | None = None,
) -> dict[str, Any]:
    """Find lattice matches/terminations without building every atomic model."""

    return copy.deepcopy(_prepare(input_path, config, cancelled)[1])


def build_frames(
    input_path: Path, config: Mapping[str, Any], *,
    progress: Callable[[int, int], None] | None = None,
    cancelled: Callable[[], bool] | None = None,
) -> dict[str, Any]:
    """Build one selected ZSL/termination pair and in-plane shift frames."""

    import numpy as np
    from pymatgen.core import Structure

    builder, search = _prepare(input_path, config, cancelled)
    index = config["match_index"]
    if not 0 <= index < len(builder.zsl_matches):
        raise ValueError("Selected match is outside the search results; search again")
    termination = tuple(config["termination"] or builder.terminations[0])
    if termination not in builder.terminations:
        raise ValueError("Selected termination is no longer available; search again")
    selected = copy.copy(builder)
    selected.zsl_matches = [builder.zsl_matches[index]]
    _check_cancelled(cancelled)
    interface = next(selected.get_interfaces(
        termination=termination, gap=config["gap"], vacuum_over_film=config["vacuum_over_film"],
        film_thickness=config["film_thickness"], substrate_thickness=config["substrate_thickness"],
        in_layers=config["in_layers"],
    ))
    _check_cancelled(cancelled)
    # Keep a/b in the interface plane and c along the slab stacking direction.
    # Plane rendering selects the normal using per-frame metadata.
    lattice = interface.lattice.matrix
    substrate_z = interface.frac_coords[interface.substrate_indices, 2]
    film_z = interface.frac_coords[interface.film_indices, 2]
    plane = float((np.max(substrate_z) + np.min(film_z)) / 2)
    base_coords = interface.frac_coords.copy()
    n1, n2 = shift_grid(config)
    total = n1 * n2
    if progress:
        progress(0, total)
    frames, metadata = [], []
    system = "-".join(structure.composition.reduced_formula for structure in (
        builder.substrate_structure, builder.film_structure,
    ))
    record = search["matches"][index]
    for i in range(n1):
        for j in range(n2):
            _check_cancelled(cancelled)
            coefficients = (
                [config["offset_a"], config["offset_b"]]
                if config["shift_mode"] == "custom" else [i / n1, j / n2]
            )
            offset = (np.asarray(coefficients) % 1).tolist()
            # Wrap only the lateral coordinates. Interface.in_plane_offset
            # also wraps z, which could split a slab crossing the cell origin.
            coords = base_coords.copy()
            coords[interface.film_indices, :2] += offset
            coords[:, :2] %= 1
            frames.append(Structure(
                lattice=lattice, species=interface.species, coords=coords,
                site_properties=interface.site_properties,
            ))
            metadata.append({
                "frame": len(frames) - 1, "engine": "pymatgen", "match_index": index,
                "shift_mode": config["shift_mode"],
                **frame_identity(config, system, coefficients, record),
                "shift_vectors_A": interface.lattice.matrix[:2].tolist(),
                "requested_shift_A": (np.asarray(coefficients) @ interface.lattice.matrix[:2]).tolist(),
                "offset_i": i, "offset_j": j, "offset_grid": [n1, n2],
                "in_plane_offset": offset, "termination": list(termination),
                "interface_plane_fractions": [plane],
                "interface_normal_axis": 2,
            })
            if progress:
                progress(len(frames), total)
    report = {
        **copy.deepcopy(search["matches"][index]), "engine": "pymatgen",
        "match_index": index, "match_count": len(search["matches"]),
        "search_truncated": search["truncated"], "termination": list(termination),
        "substrate_miller": search["substrate_miller"], "film_miller": search["film_miller"],
        "interface_area_A2": float(np.linalg.norm(np.cross(*interface.lattice.matrix[:2]))),
        "atom_count": len(interface), "frame_count": total,
        "shift_mode": config["shift_mode"],
        "system": metadata[0]["system"], "shift_basis": "supercell",
        "config": copy.deepcopy(dict(config)),
    }
    return {"frames": frames, "metadata": metadata, "report": report}

"""Bounded approximate-CSL search and atomic heterointerface construction.

InterfaceMaster supplies primitive-plane conversion, native LLL reduction and
DSC arithmetic. Exact planar integer arithmetic and the candidate workflow live
in this repository; no upstream code or working-directory files are modified.
"""

from __future__ import annotations

import copy
import math
from bisect import insort
from collections import OrderedDict
from itertools import product
from pathlib import Path
from threading import Lock
from typing import Any, Callable, Mapping

import numpy as np

from core.interface.commensurability import planar_lattices, rational_periods, reduce_planar
from core.interface.interface_backend import active_shift_basis, frame_identity, shift_grid
from core.interface.pymatgen_engine import _check_cancelled, _miller, _paths, _signature


_CACHE: OrderedDict[str, tuple[dict[str, Any], dict[str, Any]]] = OrderedDict()
_LOCK = Lock()


def _surface_frame(rows: np.ndarray) -> np.ndarray:
    a, b = rows[:2]
    x = a / np.linalg.norm(a)
    z = np.cross(a, b)
    z /= np.linalg.norm(z)
    return np.array([x, np.cross(z, x), z])


def _generator(structure: Any, miller: tuple[int, ...], thickness: float, layers: bool) -> Any:
    from pymatgen.core.surface import SlabGenerator

    # Bound the oriented cell before SlabGenerator allocates its atoms.
    divisor = math.gcd(*miller)
    reduced = [int(x // divisor) for x in miller]
    eye = np.eye(3, dtype=int)
    scale = [eye[i] for i, value in enumerate(reduced) if value == 0]
    nonzero = [i for i, value in enumerate(reduced) if value]
    multiple = math.lcm(*(reduced[i] for i in nonzero))
    for i in nonzero:
        for j in nonzero:
            if i >= j or len(scale) >= 2:
                continue
            row = [0, 0, 0]
            row[i], row[j] = -multiple // reduced[i], multiple // reduced[j]
            row_gcd = math.gcd(*row)
            scale.append(np.array(row, dtype=int) // row_gcd)
    normal = np.linalg.solve(structure.lattice.matrix, reduced)
    c_index = max(nonzero, key=lambda i: abs(np.dot(normal, structure.lattice.matrix[i])) / structure.lattice.abc[i])
    scale.append(eye[c_index])
    ouc_count = len(structure) * abs(int(round(np.linalg.det(scale))))
    if ouc_count > 200000:
        raise ValueError("The oriented surface cell exceeds 200000 atoms; choose lower plane indices or a smaller bulk input")
    return SlabGenerator(
        structure, miller, min_slab_size=thickness, min_vacuum_size=3,
        in_unit_planes=layers, primitive=False, reorient_lattice=False,
        center_slab=True,
    )


def _surfaces(input_path: Path, config: Mapping[str, Any], cancelled: Any) -> dict[str, Any]:
    from interfacemaster.cellcalc import get_primitive_hkl, get_pri_vec_inplane
    from pymatgen.core import Structure
    from pymatgen.core.interface import label_termination

    surfaces = {}
    for side, path in zip(("substrate", "film"), _paths(input_path, config)):
        _check_cancelled(cancelled)
        source = Structure.from_file(path)
        if not source.is_ordered:
            raise ValueError(f"{side}: resolve partial occupancies before constructing an atomic interface")
        input_miller = _miller(source, config, side)
        primitive = source.get_primitive_structure(tolerance=0.01, use_site_props=True)
        miller = tuple(int(x) for x in get_primitive_hkl(
            input_miller, source.lattice.matrix.T, primitive.lattice.matrix.T,
        ))
        generator = _generator(primitive, miller, 1, True)
        if len(generator.oriented_unit_cell) > 5000:
            raise ValueError("Termination enumeration exceeds 5000 oriented-cell atoms; choose lower plane indices or a smaller input cell")
        plane_vectors = get_pri_vec_inplane(miller, primitive.lattice.matrix.T)
        normal = np.linalg.solve(primitive.lattice.matrix, miller)
        if np.dot(np.cross(*plane_vectors.T), normal) < 0:
            plane_vectors[:, 1] *= -1
        rows = np.vstack((plane_vectors.T, normal))
        frame = _surface_frame(rows)
        basis = (rows @ frame.T)[:2, :2].T
        slabs = generator.get_slabs(
            ftol=config[f"{side}_termination_ftol"],
            filter_out_sym_slabs=config["filter_out_sym_slabs"],
        )
        terminations = {}
        for i, slab in enumerate(slabs, 1):
            _check_cancelled(cancelled)
            label_slab = slab.get_orthogonal_c_slab()
            label_sign = 1 if np.dot(label_slab.lattice.matrix[2], frame[2]) > 0 else -1
            if side == "film":
                label_sign *= -1
            # Identify top/bottom along the requested signed normal. Any
            # reflection here is for labeling only, never for model atoms.
            label_coords = label_slab.frac_coords.copy()
            label_coords[:, 2] = 0.5 + label_sign * (label_coords[:, 2] - 0.5)
            label_slab = Structure(label_slab.lattice, label_slab.species, label_coords)
            label = label_termination(label_slab, config[f"{side}_termination_ftol"], i if config["label_index"] else None)
            # Retain distinct cuts even if their chemical labels happen to agree.
            label = f"{label} [cut {i}]"
            terminations[label] = float(slab.shift)
        if not terminations:
            raise ValueError(f"{side}: no terminations found for the selected plane")
        surfaces[side] = {"structure": primitive, "miller": miller,
                          "input_miller": input_miller, "basis": basis,
                          "frame": frame, "plane_vectors": plane_vectors,
                          "terminations": terminations}
    return surfaces


def _integer_matrix(value: np.ndarray) -> np.ndarray:
    rounded = np.rint(value).astype(np.int64)
    if not np.allclose(value, rounded, atol=1e-6, rtol=0):
        raise ArithmeticError("CSL vectors are not integral in both surface lattices")
    return rounded


def _verify_native_dscl(record: dict[str, Any]) -> str:
    """Cross-check the exact planar sum with InterfaceMaster's DSC solver.

    A common normal is appended to the planar problem, rather than assuming
    that a projection of a three-dimensional DSC is the planar DSCL.
    """

    from interfacemaster.cellcalc import DSCcalc

    rational = np.asarray(record["rational_numerator"]) / record["rational_denominator"]
    embedded = np.eye(3)
    embedded[:2, :2] = rational
    try:
        calculator = DSCcalc()
        calculator.parse_int_U(np.eye(3), embedded, record["rational_denominator"], tol=1e-7)
        calculator.compute_DSC(to_LLL=False, tol=1e-7)
        native = np.asarray(calculator.DSC)
        if not np.allclose(native[2, :2], 0, atol=1e-7):
            raise ArithmeticError("native DSC mixed the auxiliary normal into the plane")
        basis = native[:2, :2]
        _integer_matrix(np.linalg.solve(basis, np.eye(2)))
        _integer_matrix(np.linalg.solve(basis, rational))
        substrate = np.asarray(record["substrate_primitive_vectors_A"]).T
        area = abs(float(np.linalg.det(substrate @ basis)))
        if not math.isclose(area, record["dscl_area_A2"], rel_tol=1e-6):
            raise ArithmeticError("native DSC area differs from the exact planar result")
    except (ValueError, RuntimeError, ArithmeticError, np.linalg.LinAlgError) as exc:
        # Preserve a usable exact result if native floating arithmetic fails.
        # Keep the reason in the report so this is never an invisible fallback.
        return f"Exact planar HNF; native DSC cross-check unavailable: {exc}"
    return "Exact planar HNF; independently confirmed by InterfaceMaster DSCcalc"


def _candidate(substrate: np.ndarray, film: np.ndarray, angle: float,
               numerator: np.ndarray, denominator: int, config: Mapping[str, Any]) -> dict[str, Any] | None:
    from pymatgen.analysis.elasticity.strain import Deformation

    radians = math.radians(angle)
    rotation = np.array([[math.cos(radians), -math.sin(radians)],
                         [math.sin(radians), math.cos(radians)]])
    target = substrate @ (numerator / denominator)
    deformation = target @ np.linalg.inv(rotation @ film)
    left, stretches, right = np.linalg.svd(deformation)
    principal = stretches - 1
    if np.max(np.abs(principal)) > config["im_max_strain"] + 1e-10:
        return None
    if abs(np.linalg.det(deformation) - 1) > config["im_max_area_change"] + 1e-10:
        return None
    residual_rotation = left @ right
    correction = math.degrees(math.atan2(residual_rotation[1, 0], residual_rotation[0, 0]))
    actual_angle = angle + correction
    if (abs(correction) > config["im_rotation_tol"] + 1e-8
            or actual_angle < config["im_angle_start"] - 1e-7
            or actual_angle > config["im_angle_stop"] + 1e-7):
        return None
    fractional_csl, fractional_dscl = planar_lattices(numerator, denominator)
    csl = reduce_planar(substrate @ fractional_csl)
    sub_matrix = _integer_matrix(np.linalg.solve(substrate, csl))
    film_matrix = _integer_matrix(np.linalg.solve(target, csl))
    sub_index = int(round(abs(np.linalg.det(sub_matrix))))
    film_index = int(round(abs(np.linalg.det(film_matrix))))
    sigma = math.sqrt(sub_index * film_index)
    sub_area = abs(float(np.linalg.det(csl)))
    raw_film = rotation @ film @ film_matrix
    film_area = abs(float(np.linalg.det(raw_film)))
    if (max(sub_index, film_index) > config["im_max_index"]
            or sigma > config["im_max_sigma"] + 1e-8
            or max(sub_area, film_area) > config["max_area"] + 1e-8):
        return None
    dscl = reduce_planar(substrate @ fractional_dscl)
    dscl_area = abs(float(np.linalg.det(dscl)))
    if not math.isclose(math.sqrt(sub_area / dscl_area), sigma, rel_tol=1e-7):
        raise ArithmeticError("CSL/DSCL area identity failed")
    total_deformation = np.eye(3)
    total_deformation[:2, :2] = target @ np.linalg.inv(film)
    strain = Deformation(total_deformation).green_lagrange_strain

    def internal_angle(vectors: np.ndarray) -> float:
        return math.atan2(abs(float(np.linalg.det(vectors))), float(np.dot(*vectors.T)))

    film_angle, sub_angle = internal_angle(raw_film), internal_angle(csl)
    return {
        "index": 0, "twist_angle_deg": actual_angle, "sample_angle_deg": angle,
        "rotation_correction_deg": correction,
        "substrate_index": sub_index, "film_index": film_index, "sigma": sigma,
        "substrate_area_A2": sub_area, "film_area_A2": film_area,
        "dscl_area_A2": dscl_area,
        "substrate_supercell_matrix": sub_matrix.T.tolist(),
        "film_supercell_matrix": film_matrix.T.tolist(),
        "substrate_vectors_A": np.column_stack((csl.T, np.zeros(2))).tolist(),
        "film_vectors_A": np.column_stack((raw_film.T, np.zeros(2))).tolist(),
        "substrate_primitive_vectors_A": substrate.T.tolist(),
        "film_approximant_vectors_A": target.T.tolist(),
        "dscl_vectors_A": dscl.T.tolist(),
        "principal_strains": principal.tolist(),
        "von_mises_strain": float(strain.von_mises_strain), "strain": np.asarray(strain).tolist(),
        "deformation": total_deformation.tolist(),
        "length_mismatch": (np.linalg.norm(csl, axis=0) / np.linalg.norm(raw_film, axis=0) - 1).tolist(),
        "angle_mismatch": sub_angle / film_angle - 1,
        "angle_difference_deg": math.degrees(sub_angle - film_angle),
        "rational_numerator": numerator.tolist(), "rational_denominator": denominator,
    }


def _rank(record: dict[str, Any], key: str) -> tuple[float, ...]:
    strain = max(abs(x) for x in record["principal_strains"])
    if key == "strain":
        return strain, record["sigma"], record["substrate_area_A2"]
    if key == "area":
        return record["substrate_area_A2"], strain, record["sigma"]
    return record["sigma"], strain, record["substrate_area_A2"]


def _prepare(input_path: Path, config: Mapping[str, Any], cancelled: Any) -> tuple[dict[str, Any], dict[str, Any]]:
    _check_cancelled(cancelled)
    signature = _signature(input_path, config)
    with _LOCK:
        if signature in _CACHE:
            _CACHE.move_to_end(signature)
            return _CACHE[signature]
    surfaces = _surfaces(input_path, config, cancelled)
    substrate, film = (surfaces[side]["basis"] for side in ("substrate", "film"))
    count = math.ceil((config["im_angle_stop"] - config["im_angle_start"]) / config["im_angle_step"])
    retained, seen, accepted = [], set(), 0
    for i in range(count + 1):
        angle = min(config["im_angle_start"] + i * config["im_angle_step"], config["im_angle_stop"])
        radians = math.radians(angle)
        rotation = np.array([[math.cos(radians), -math.sin(radians)], [math.sin(radians), math.cos(radians)]])
        transform = np.linalg.solve(substrate, rotation @ film)
        for denominator in range(1, config["im_max_denominator"] + 1):
            _check_cancelled(cancelled)
            numerator = np.rint(transform * denominator).astype(np.int64)
            divisor = math.gcd(denominator, *(int(x) for x in numerator.ravel()))
            numerator, q = numerator // divisor, denominator // divisor
            key = (*numerator.ravel().tolist(), q)
            determinant = int(numerator[0, 0]) * int(numerator[1, 1]) - int(numerator[0, 1]) * int(numerator[1, 0])
            if key in seen or determinant <= 0:
                continue
            record = _candidate(substrate, film, angle, numerator, q, config)
            if record is None:
                continue
            seen.add(key)
            accepted += 1
            entry = (_rank(record, config["im_rank_by"]), accepted, record)
            if len(retained) < config["max_matches"] or entry[:2] < retained[-1][:2]:
                insort(retained, entry)
                if len(retained) > config["max_matches"]:
                    retained.pop()
    matches = [row for _rank_key, _order, row in retained]
    if not matches:
        raise ValueError("No approximate CSL found within the rotation, denominator, strain, index and area limits. Increase bounds or relax strain limits.")
    for i, record in enumerate(matches):
        record["index"] = i
    if len(surfaces["film"]["terminations"]) * len(surfaces["substrate"]["terminations"]) > 10000:
        raise ValueError("More than 10000 termination pairs; increase plane grouping tolerances or enable symmetry filtering")
    terminations = list(product(surfaces["film"]["terminations"], surfaces["substrate"]["terminations"]))
    result = {
        "engine": "interface-master", "matches": matches,
        "terminations": [list(pair) for pair in terminations],
        "truncated": accepted > len(matches), "total_matches": accepted,
        "trial_count": (count + 1) * config["im_max_denominator"],
        "substrate_miller": list(surfaces["substrate"]["input_miller"]),
        "film_miller": list(surfaces["film"]["input_miller"]),
        "primitive_millers": {side: list(surfaces[side]["miller"]) for side in surfaces},
        "rational_periods": {
            axis: rational_periods(float(np.linalg.norm(film[:, j])), float(np.linalg.norm(substrate[:, j])), config["im_max_denominator"])
            for j, axis in enumerate(("a", "b"))
        },
    }
    _check_cancelled(cancelled)
    with _LOCK:
        _CACHE[signature] = (surfaces, result)
        _CACHE.move_to_end(signature)
        while len(_CACHE) > 2:
            _CACHE.popitem(last=False)
    return surfaces, result


def search_matches(input_path: Path, config: Mapping[str, Any], *, cancelled: Any = None) -> dict[str, Any]:
    """Retain the best candidates over the whole requested finite search."""

    return copy.deepcopy(_prepare(input_path, config, cancelled)[1])


def _matched_slab(surface: dict[str, Any], label: str, side: str,
                  record: dict[str, Any], config: Mapping[str, Any]) -> Any:
    from pymatgen.core import Structure

    generator = _generator(surface["structure"], surface["miller"], config[f"{side}_thickness"], config["in_layers"])
    ouc = generator.oriented_unit_cell
    normal = surface["frame"][2]
    ouc_height = abs(float(np.dot(ouc.lattice.matrix[2], normal)))
    requested_height = (config[f"{side}_thickness"] * surface["structure"].lattice.d_hkl(surface["miller"])
                        if config["in_layers"] else config[f"{side}_thickness"])
    if len(ouc) * math.ceil(requested_height / ouc_height) > 200000:
        raise ValueError("The oriented slab intermediate exceeds 200000 atoms; reduce thickness or plane indices")
    slab = generator.get_slab(
        shift=surface["terminations"][label],
    ).get_orthogonal_c_slab()
    # SlabGenerator's oriented cell can contain several primitive surface
    # cells. Replace that basis explicitly, then remove only exact replicas.
    # This keeps CSL indices independent of an incidental slab supercell.
    primitive_rows = np.vstack((surface["plane_vectors"].T,
                                normal * abs(float(np.dot(slab.lattice.matrix[2], normal)))))
    ratio = slab.volume / abs(float(np.linalg.det(primitive_rows)))
    multiplicity = int(round(ratio))
    if multiplicity < 1 or not math.isclose(ratio, multiplicity, rel_tol=1e-7):
        raise ArithmeticError("The oriented slab is not integral in the primitive surface basis")
    frac = slab.cart_coords @ np.linalg.inv(primitive_rows)
    frac[:, :2] %= 1
    keys = np.round(frac, 8)
    keys[:, :2] %= 1
    unique = {}
    for i, site in enumerate(slab):
        unique.setdefault((site.species_string, *keys[i]), i)
    indices = list(unique.values())
    if len(indices) * multiplicity != len(slab):
        raise ArithmeticError("Atomic sites do not reproduce the primitive surface translations")
    slab = Structure(primitive_rows, [slab[i].specie for i in indices], frac[indices],
                     site_properties={key: [values[i] for i in indices] for key, values in slab.site_properties.items()})
    # Apply integer repetitions before setting the common lateral lattice.
    transform = np.eye(3, dtype=int)
    transform[:2, :2] = record[f"{side}_supercell_matrix"]
    expected = len(slab) * abs(int(round(np.linalg.det(transform))))
    if expected > config["im_max_atoms"]:
        raise ValueError(f"{side}: {expected} atoms exceeds the configured atom limit")
    slab.make_supercell(transform)
    frame = surface["frame"]
    coords = slab.cart_coords @ frame.T
    rows = slab.lattice.matrix @ frame.T
    if rows[2, 2] <= 0:
        raise ValueError("Slab normal points away from the selected interface surface")
    # Each slab keeps its own in-plane fractional coordinates; replacing film
    # rows therefore applies exactly the deformation accepted by the search.
    frac = coords @ np.linalg.inv(rows)
    common = np.asarray(record["substrate_vectors_A"])
    lattice = np.vstack((common, [0, 0, rows[2, 2]]))
    return Structure(lattice, slab.species, frac, site_properties=slab.site_properties)


def build_frames(input_path: Path, config: Mapping[str, Any], *,
                 progress: Callable[[int, int], None] | None = None,
                 cancelled: Callable[[], bool] | None = None) -> dict[str, Any]:
    """Build explicit terminations and mutually exclusive custom/sample shifts."""

    from pymatgen.core import Structure

    from core.interface.registry_analysis import contact_geometry

    surfaces, search = _prepare(input_path, config, cancelled)
    index = config["match_index"]
    if index >= len(search["matches"]):
        raise ValueError("The selected CSL candidate is unavailable; search again")
    record = search["matches"][index]
    dscl_method = _verify_native_dscl(record)
    termination = config["termination"] or search["terminations"][0]
    if termination not in search["terminations"]:
        raise ValueError("The selected terminations are unavailable; search again")
    _check_cancelled(cancelled)
    substrate = _matched_slab(surfaces["substrate"], termination[1], "substrate", record, config)
    film = _matched_slab(surfaces["film"], termination[0], "film", record, config)
    if len(substrate) + len(film) > config["im_max_atoms"]:
        raise ValueError("The combined model exceeds im_max_atoms; reduce thickness or choose a smaller CSL")
    common = np.asarray(record["substrate_vectors_A"])[:, :2]
    sub_coords, film_coords = substrate.cart_coords.copy(), film.cart_coords.copy()
    sub_coords[:, 2] -= sub_coords[:, 2].min()
    # Film bottom meets substrate top; retain the film's handedness.
    film_coords[:, 2] -= film_coords[:, 2].min()
    film_coords[:, 2] += sub_coords[:, 2].max() + config["gap"]
    plane_z = (sub_coords[:, 2].max() + film_coords[:, 2].min()) / 2
    height = film_coords[:, 2].max() + config["vacuum_over_film"]
    if height <= 1e-8:
        raise ValueError("The interface has zero height; increase thickness, gap or vacuum")
    z_shift = config["vacuum_over_film"] / 2
    sub_coords[:, 2] += z_shift
    film_coords[:, 2] += z_shift
    lattice = np.array([[*common[0], 0], [*common[1], 0], [0, 0, height]])
    coords = np.vstack((sub_coords, film_coords)) @ np.linalg.inv(lattice)
    species = substrate.species + film.species
    properties = {
        key: substrate.site_properties[key] + film.site_properties[key]
        for key in set(substrate.site_properties) & set(film.site_properties)
    }
    properties["interface_label"] = ["substrate"] * len(substrate) + ["film"] * len(film)
    n1, n2 = shift_grid(config)
    total = n1 * n2
    if len(species) * total > 2000000:
        raise ValueError("The sampled frames exceed two million total atoms; reduce the grid or model size")
    basis_name = active_shift_basis(config)
    sample_basis = np.asarray(record["dscl_vectors_A"]) if basis_name == "dscl" else common
    sample_in_common = sample_basis @ np.linalg.inv(common) if basis_name == "dscl" else np.eye(2)
    dscl_fractional = np.column_stack((
        np.asarray(record["dscl_vectors_A"]) @ np.linalg.inv(common), np.zeros(2),
    )).tolist()
    frames, metadata = [], []
    system = "-".join(surfaces[side]["structure"].composition.reduced_formula for side in ("substrate", "film"))
    if progress:
        progress(0, total)
    for i in range(n1):
        for j in range(n2):
            _check_cancelled(cancelled)
            coefficients = (np.array([config["offset_a"], config["offset_b"]]) if config["shift_mode"] == "custom"
                            else np.array([i / n1, j / n2]))
            offset = coefficients @ sample_in_common
            shifted = coords.copy()
            shifted[len(substrate):, :2] += offset
            shifted[:, :2] %= 1
            geometry = contact_geometry(
                shifted @ lattice, len(substrate), common, config, cancelled,
            )
            frames.append(Structure(lattice, species, shifted, site_properties=properties))
            metadata.append({
                "engine": "interface-master", "frame": len(frames) - 1,
                "match_index": index, "shift_mode": config["shift_mode"],
                "twist_angle_deg": record["twist_angle_deg"], "sigma": record["sigma"],
                "substrate_index": record["substrate_index"], "film_index": record["film_index"],
                "dscl_area_A2": record["dscl_area_A2"],
                "dscl_vectors_fractional": dscl_fractional,
                **frame_identity(config, system, coefficients.tolist(), record),
                "shift_vectors_A": np.column_stack((sample_basis, np.zeros(2))).tolist(),
                "requested_shift_A": [*(coefficients @ sample_basis).tolist(), 0.0],
                "offset_i": i, "offset_j": j, "offset_grid": [n1, n2],
                "in_plane_offset": (offset % 1).tolist(), "termination": list(termination),
                "interface_plane_fractions": [(plane_z + z_shift) / height],
                "interface_normal_axis": 2, **geometry,
            })
            if progress:
                progress(len(frames), total)
    best = min(metadata, key=lambda row: (-row["registry_fraction"], row["mean_lateral_distance_A"]))
    report = {
        **copy.deepcopy(record), "engine": "interface-master", "match_index": index,
        "match_count": len(search["matches"]), "search_truncated": search["truncated"],
        "termination": list(termination), "substrate_miller": search["substrate_miller"],
        "film_miller": search["film_miller"], "interface_area_A2": record["substrate_area_A2"],
        "atom_count": len(species), "frame_count": total, "shift_mode": config["shift_mode"],
        "shift_basis": best["shift_basis"], "registry_best_frame": best["frame"],
        "system": best["system"],
        "registry_best": best, "clashing_frames": sum(row["close_contact_count"] > 0 for row in metadata),
        "shift_frames": copy.deepcopy(metadata),
        "rational_periods": copy.deepcopy(search["rational_periods"]),
        "dscl_method": dscl_method,
        "config": copy.deepcopy(dict(config)),
    }
    return {"frames": frames, "metadata": metadata, "report": report}

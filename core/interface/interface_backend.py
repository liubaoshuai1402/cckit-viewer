# CCkit 功能说明：
# 使用 InterfaceMaster 匹配两种材料的接触面并生成 CNID 网格结构。
"""Interface Builder backend for two structures and a two-dimensional CSL.

The substrate remains the reference lattice. InterfaceMaster aligns the film,
searches the requested twist interval, and applies its accepted transformation
to the film. The first admissible match is used for all CNID sampling frames.
This module is a desktop feature adapter; it does not replace model_interface.py.
"""

from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import json
import math
import os
import tempfile
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Mapping


@lru_cache(maxsize=1)
def _gb() -> Any:
    """Share OmniGB's validated inputs and InterfaceMaster dependency imports."""

    path = Path(__file__).resolve().parents[1] / "omnigb" / "omnigb_backend.py"
    spec = importlib.util.spec_from_file_location("cckit_omnigb_backend", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"could not load OmniGB backend: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def default_config() -> dict[str, Any]:
    """Return an independent editable configuration for the desktop card."""

    return {
        "mode": "build",
        "film_structure": "",
        "substrate_plane": [0, 0, 1],
        "substrate_plane_notation": "three",
        "film_plane": [0, 0, 1],
        "film_plane_notation": "three",
        "align_directions": False,
        "substrate_direction": [1, 0, 0],
        "substrate_direction_notation": "three",
        "film_direction": [1, 0, 0],
        "film_direction_notation": "three",
        "angle_start_deg": 0.0,
        "angle_stop_deg": 180.0,
        "angle_step_deg": 1.0,
        "max_sigma_substrate": 100,
        "max_sigma_film": 100,
        "csl_mode": "near",
        "tolerances": {"du": 0.02, "area": 0.02, "deformation": 0.02, "integer": 1e-3},
        "substrate_normal_repeat": 1,
        "film_normal_repeat": 1,
        "normal_ortho": False,
        "superlattice_lim": 50,
        "dp1": 0.0,
        "dp2": 0.0,
        "normal_distance": 0.0,
        "vacuum": 0.0,
        "KTI": False,
        "cnid_grid": [4, 4],
    }


def parameters(config: Mapping[str, Any]) -> dict[str, Any]:
    """Validate desktop values without running a CSL search or reading files."""

    gb = _gb()
    values = default_config()
    values.update(config)
    if values["mode"] != "build":
        raise ValueError("Interface Builder supports build mode only")
    values["tolerances"] = gb._tolerances(values)
    values["csl_mode"], _ = gb._csl_mode(values)
    for prefix in ("substrate", "film"):
        for kind in ("plane", "direction"):
            key = f"{prefix}_{kind}"
            notation = gb.index_notation(values[key], config.get(f"{key}_notation"))
            values[f"{key}_notation"] = notation
            if kind == "plane" or values["align_directions"]:
                values[key] = gb.index_values(values[key], kind, notation)
        for key in (f"max_sigma_{prefix}", f"{prefix}_normal_repeat"):
            values[key] = gb._positive_int(values[key], key)
    for key in ("angle_start_deg", "angle_stop_deg", "angle_step_deg", "dp1", "dp2"):
        values[key] = gb._finite_float(values[key], key)
    if values["angle_stop_deg"] < values["angle_start_deg"]:
        raise ValueError("angle_stop_deg must be greater than or equal to angle_start_deg")
    if values["angle_step_deg"] <= 0:
        raise ValueError("angle_step_deg must be positive")
    for key in ("normal_distance", "vacuum"):
        values[key] = gb._nonnegative_float(values[key], key)
    for key in ("align_directions", "normal_ortho", "KTI"):
        if not isinstance(values[key], bool):
            raise ValueError(f"{key} must be a JSON boolean")
    values["superlattice_lim"] = gb._positive_int(values["superlattice_lim"], "superlattice_lim")
    values["cnid_grid"] = list(gb.cnid_grid(values))
    return values


def build_interface_frames(
    input_path: Path,
    config: Mapping[str, Any],
    *,
    progress: Callable[[int, int], None] | None = None,
    cancelled: Callable[[], bool] | None = None,
) -> dict[str, Any]:
    """Match specified surfaces and sample the resulting CNID with native APIs."""

    gb = _gb()
    values = parameters(config)
    np, Structure, primitive_hkl, _, convert_direction, core = gb._backend()
    from interfacemaster.interface_generator import get_disorientation, get_height
    from interfacemaster.cellcalc import get_normal_from_MI

    substrate_path = input_path.expanduser().resolve()
    if not substrate_path.is_file():
        raise FileNotFoundError(f"substrate structure was not found: {substrate_path}")
    film_value = values.get("film_structure")
    if not isinstance(film_value, str) or not film_value.strip():
        raise ValueError("select a film structure before building the interface")
    film_path = Path(film_value).expanduser()
    if not film_path.is_absolute():
        film_path = substrate_path.parent / film_path
    film_path = film_path.resolve()
    if not film_path.is_file():
        raise FileNotFoundError(f"film structure was not found: {film_path}")
    substrate = Structure.from_file(substrate_path)
    film = Structure.from_file(film_path)
    interface = core(substrate, film, verbose=False)
    planes = []
    directions = []
    for prefix, original, primitive in (
        ("substrate", interface.conv_lattice_1, interface.lattice_1),
        ("film", interface.conv_lattice_2, interface.lattice_2),
    ):
        plane = gb.three_indices(
            values[f"{prefix}_plane"], "plane", values[f"{prefix}_plane_notation"], original
        )
        plane = primitive_hkl(plane, original, primitive, values["tolerances"]["integer"])
        planes.append(plane)
        if values["align_directions"]:
            direction = gb.three_indices(
                values[f"{prefix}_direction"], "direction",
                values[f"{prefix}_direction_notation"], original,
            )
            direction = convert_direction(original, primitive, np.asarray(direction))
            normal = get_normal_from_MI(primitive, plane)
            cart = primitive @ direction
            if abs(np.dot(normal, cart)) > 1e-7 * np.linalg.norm(normal) * np.linalg.norm(cart):
                raise ValueError(f"{prefix} reference direction must lie in its contact plane")
            directions.append(direction)
    tolerances = values["tolerances"]
    interface.parse_limit(
        tolerances["du"], tolerances["area"], values["max_sigma_substrate"],
        values["max_sigma_film"], tolerances["deformation"],
    )
    orientation_args = {}
    if values["align_directions"]:
        orientation_args = {
            "pre_dt": True,
            "pre_R": get_disorientation(
                interface.lattice_1, interface.lattice_2,
                directions[0], planes[0], directions[1], planes[1],
            ),
        }
    n1, n2 = values["cnid_grid"]
    total = n1 * n2
    if progress is not None:
        progress(0, total)
    frames, metadata = [], []
    with tempfile.TemporaryDirectory(prefix="cckit_heterointerface_") as work_dir:
        old_cwd = Path.cwd()
        try:
            os.chdir(work_dir)
            if cancelled is not None and cancelled():
                raise RuntimeError("interface construction was cancelled")
            # A fresh object has a placeholder CNID. Clear it so a failed
            # native search cannot be mistaken for a successful match.
            interface.CNID = None
            # Search one native angle at a time so superseded GUI requests can
            # stop between angles. This also includes the requested stop angle.
            for angle in gb._scan_angles(
                values["angle_start_deg"], values["angle_stop_deg"], values["angle_step_deg"]
            ):
                if cancelled is not None and cancelled():
                    raise RuntimeError("interface construction was cancelled")
                with contextlib.redirect_stdout(io.StringIO()):
                    interface.search_one_position_2D(
                        planes[0], planes[1],
                        theta_range=values["angle_step_deg"],
                        dtheta=values["angle_step_deg"], start=angle,
                        integer_tol=tolerances["integer"],
                        exact=values["csl_mode"] == "exact", **orientation_args,
                    )
                if interface.CNID is not None:
                    break
            if interface.CNID is None:
                raise RuntimeError(
                    "No 2D CSL match found. Adjust the twist interval, maximum "
                    "surface indices or mismatch tolerances."
                )
            interface.compute_bicrystal_two_D(
                hkl_1=planes[0], hkl_2=planes[1],
                lim=values["superlattice_lim"], normal_ortho=values["normal_ortho"],
                tol_ortho=tolerances["integer"], tol_integer=tolerances["integer"],
            )
            vectors = np.asarray(interface.CNID, dtype=float)
            if vectors.shape != (3, 2) or not np.isfinite(vectors).all():
                raise RuntimeError("InterfaceMaster returned an invalid CNID")
            if np.linalg.norm(np.cross(vectors[:, 0], vectors[:, 1])) < 1e-12:
                raise RuntimeError("InterfaceMaster returned a degenerate CNID")
            for i in range(n1):
                for j in range(n2):
                    if cancelled is not None and cancelled():
                        raise RuntimeError("interface construction was cancelled")
                    alpha, beta = i / n1, j / n2
                    translation = vectors @ np.asarray([alpha, beta])
                    interface.get_bicrystal(
                        two_D=True, dydz=translation,
                        xyz_1=[values["substrate_normal_repeat"], 1, 1],
                        xyz_2=[values["film_normal_repeat"], 1, 1],
                        dp1=values["dp1"], dp2=values["dp2"],
                        dx=values["normal_distance"], vx=values["vacuum"],
                        KTI=values["KTI"], filename=str(Path(work_dir) / "POSCAR_frame"),
                        filetype="VASP",
                    )
                    frames.append(interface.bicrystal_structure.copy())
                    height = get_height(interface.lattice_bi)
                    contact = (get_height(interface.slab_lattice_1) + values["normal_distance"] / 2) / height
                    fractions = [contact]
                    if values["vacuum"] == 0:
                        fractions.append(1 - values["normal_distance"] / (2 * height))
                    metadata.append({
                        "frame": len(frames) - 1, "cnid_i": i, "cnid_j": j,
                        "cnid_n1": n1, "cnid_n2": n2,
                        "rbt_alpha": alpha, "rbt_beta": beta,
                        "rbt_translation_A": translation.tolist(),
                        "cnid_u1_A": vectors[:, 0].tolist(), "cnid_u2_A": vectors[:, 1].tolist(),
                        "termination_dp1_A": values["dp1"], "termination_dp2_A": values["dp2"],
                        "interface_plane_fractions": fractions,
                    })
                    if progress is not None:
                        progress(len(frames), total)
        finally:
            os.chdir(old_cwd)
    plane_vectors = interface.lattice_1 @ interface.U1
    report = {
        "substrate_structure": str(substrate_path), "film_structure": str(film_path),
        "substrate_plane_primitive": np.asarray(planes[0]).tolist(),
        "film_plane_primitive": np.asarray(planes[1]).tolist(),
        "twist_angle_deg": float(interface.theta * 180 / np.pi),
        "substrate_area_index": int(round(abs(np.linalg.det(interface.cell_calc.U1)))),
        "film_area_index": int(round(abs(np.linalg.det(interface.cell_calc.U2)))),
        "interface_area_A2": float(np.linalg.norm(np.cross(plane_vectors[:, 0], plane_vectors[:, 1]))),
        "D": np.asarray(interface.D).tolist(),
        "deformation_max_abs": float(np.max(np.abs(interface.D - np.eye(3)))),
        "film_transform": np.asarray(interface.a2_transform).tolist(),
        "film_metric_strain_max": float(np.max(np.abs(
            np.linalg.svd(interface.a2_transform, compute_uv=False) - 1
        ))),
        "atom_count": len(frames[0]), "frame_count": len(frames),
        "csl_mode": values["csl_mode"], "selection": "first admissible 2D CSL match",
        "config": values,
    }
    return {"frames": frames, "metadata": metadata, "report": report}


def _default_output_name(mode: str, config: Mapping[str, Any]) -> str:
    """Give exported frame sequences a stable interface-specific basename."""

    return "POSCAR_Interface"


def main() -> None:
    """Allow direct developer use without registering another public command."""

    parser = argparse.ArgumentParser(description="Build a two-material interface using InterfaceMaster.")
    parser.add_argument("substrate", type=Path, help="substrate structure")
    parser.add_argument("film", type=Path, help="film structure")
    parser.add_argument("--config", type=Path, required=True, help="Interface Builder JSON settings")
    parser.add_argument("-o", "--output", type=Path, help="POSCAR output; bare names use substrate directory")
    parser.add_argument("--frame", type=int, default=0, help="zero-based CNID frame index to export")
    args = parser.parse_args()
    substrate = args.substrate.expanduser().resolve()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        parser.error("configuration must be a JSON object")
    config["film_structure"] = str(args.film.expanduser().resolve())
    output = args.output or Path("POSCAR_Interface")
    if not output.is_absolute() and output.parent == Path("."):
        output = substrate.parent / output
    output = output.expanduser().resolve()
    try:
        result = build_interface_frames(substrate, config)
        if not 0 <= args.frame < len(result["frames"]):
            raise ValueError("frame index is outside the sampled grid")
        output.parent.mkdir(parents=True, exist_ok=True)
        result["frames"][args.frame].to(fmt="poscar", filename=str(output))
        output.with_name(output.name + ".json").write_text(
            json.dumps(result["report"], indent=2), encoding="utf-8"
        )
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    print(f"Exported frame {args.frame} of {len(result['frames'])}: {output}")


if __name__ == "__main__":
    main()

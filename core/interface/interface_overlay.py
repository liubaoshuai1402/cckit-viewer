"""Render a DSCL unit cell on the built atomic interface plane."""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np


def add_dscl_overlay(plotter: Any, pv: Any, structure: Any,
                     metadata: Mapping[str, Any]) -> None:
    """Draw actual-length DSCL vectors in the current structure's coordinates.

    The normal axis selects the constant fractional coordinate of the plane.
    New frames use a/b in-plane and c normal; older frames without this marker
    retain their first-axis convention. Fractional vectors follow the frame's
    lattice and remain valid under rigid reorientation.
    """

    if metadata.get("engine") != "interface-master":
        return
    try:
        fractional = np.asarray(metadata.get("dscl_vectors_fractional", []), dtype=float)
        planes = np.asarray(metadata.get("interface_plane_fractions", []), dtype=float)
        normal_axis = int(metadata.get("interface_normal_axis", 0))
    except (TypeError, ValueError, OverflowError):
        return
    if (normal_axis not in (0, 1, 2)
            or fractional.shape != (2, 3) or planes.ndim != 1 or not len(planes)
            or not np.all(np.isfinite(fractional)) or not np.all(np.isfinite(planes))
            or not np.allclose(fractional[:, normal_axis], 0, atol=1e-10, rtol=0)):
        return
    lattice = np.asarray(structure.lattice.matrix, dtype=float)
    vectors = fractional @ lattice
    lengths = np.linalg.norm(vectors, axis=1)
    if np.any(lengths <= 1e-10):
        return
    center_fractional = np.full(3, 0.5)
    center_fractional[normal_axis] = planes[0]
    center = center_fractional @ lattice
    origin = center - 0.5 * vectors.sum(axis=0)
    vertices = np.array([origin, origin + vectors[0],
                         origin + vectors.sum(axis=0), origin + vectors[1]])
    color = "#9333ea"
    outline = pv.PolyData(vertices, lines=np.array([5, 0, 1, 2, 3, 0]))
    plotter.add_mesh(outline, color=color, line_width=5, render_lines_as_tubes=True,
                     name="interface_dscl_cell", reset_camera=False)
    for i, vector in enumerate(vectors):
        arrow = pv.Arrow(start=origin, direction=vector, scale="auto",
                         shaft_radius=0.015, tip_radius=0.045)
        plotter.add_mesh(arrow, color=color, name=f"interface_dscl_d{i + 1}",
                         reset_camera=False)
    plotter.add_point_labels(
        origin + vectors, ["DSCL d1", "DSCL d2"], text_color=color,
        font_size=12, show_points=False, shape=None, always_visible=True,
        name="interface_dscl_vectors", reset_camera=False,
    )
    plotter.add_text(
        f"DSCL on interface (actual scale)\n|d1| {lengths[0]:.3f} Å  |d2| {lengths[1]:.3f} Å",
        position="upper_right", font_size=10, color=color, name="interface_dscl_label",
    )

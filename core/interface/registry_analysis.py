"""Species-independent contact geometry with exact lateral periodic images."""

from __future__ import annotations

import math
from typing import Any, Callable, Mapping

import numpy as np
from scipy.spatial import cKDTree


def _nearest_periodic(query: np.ndarray, targets: np.ndarray, common: np.ndarray,
                      cancelled: Callable[[], bool] | None, *, all_distances: bool,
                      threshold: float = 0) -> np.ndarray:
    """Expand image shells until a singular-value bound excludes every unseen image.

    Query and target lateral fractional coordinates lie in [0, 1). An unseen
    image outside shell k has at least one fractional separation >= k, so
    its distance is at least k times the smallest lattice singular value.
    In summary mode, only the global minimum and distances below threshold
    are guaranteed; in all-distances mode every returned minimum is exact.
    """

    tree = cKDTree(targets)
    best = np.full(len(query), np.inf)
    minimum_stretch = np.linalg.svd(common, compute_uv=False).min()
    normal_separation = 0.0
    if query.shape[1] == 3:
        normal_separation = max(0.0, float(query[:, 2].min() - targets[:, 2].max()),
                                float(targets[:, 2].min() - query[:, 2].max()))
    if minimum_stretch <= 1e-10:
        raise ValueError("The matched surface cell is numerically singular")
    for shell in range(65):
        for i in range(-shell, shell + 1):
            for j in range(-shell, shell + 1):
                if max(abs(i), abs(j)) != shell:
                    continue
                if cancelled is not None and cancelled():
                    raise RuntimeError("Interface operation cancelled")
                translation = np.zeros(query.shape[1])
                translation[:2] = np.array([i, j]) @ common
                distances = tree.query(query - translation)[0]
                np.minimum(best, distances, out=best)
        bound = float(best.max()) if all_distances else max(float(best.min()), threshold)
        if math.hypot(minimum_stretch * shell, normal_separation) > bound:
            return best
    raise ValueError("Periodic contact search could not bound the oblique cell; select a less elongated CSL")


def contact_geometry(coords: np.ndarray, substrate_count: int, common: np.ndarray,
                     config: Mapping[str, Any], cancelled: Any) -> dict[str, Any]:
    """Describe registry of contact-layer atoms, without energetic ranking.

    Registry is directed: each film-layer atom finds the nearest substrate
    projection. Several film atoms may share one substrate neighbor.
    Close contacts count film atoms, rather than all cross-interface pairs.
    """

    substrate, film = coords[:substrate_count], coords[substrate_count:]
    sub_layer = substrate[substrate[:, 2] >= substrate[:, 2].max() - config["im_contact_depth"] - 1e-8]
    film_layer = film[film[:, 2] <= film[:, 2].min() + config["im_contact_depth"] + 1e-8]
    lateral = _nearest_periodic(film_layer[:, :2], sub_layer[:, :2], common, cancelled, all_distances=True)
    physical = _nearest_periodic(film, substrate, common, cancelled, all_distances=False, threshold=config["im_min_distance"])
    return {
        "contact_film_atoms": len(film_layer), "contact_substrate_atoms": len(sub_layer),
        "registry_fraction": float(np.mean(lateral <= config["im_registry_tolerance"] + 1e-8)),
        "mean_lateral_distance_A": float(np.mean(lateral)),
        "minimum_cross_distance_A": float(physical.min()),
        "close_contact_count": int(np.count_nonzero(physical < config["im_min_distance"])),
    }

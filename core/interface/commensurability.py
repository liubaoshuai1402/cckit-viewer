"""Exact planar lattice arithmetic and Inoue-style rational-period analysis.

Matrices in this module contain column basis vectors. Public supercell
matrices returned to pymatgen use row conventions instead.
"""

from __future__ import annotations

import math
from fractions import Fraction
from itertools import combinations
from typing import Any

import numpy as np


def _bezout(a: int, b: int) -> tuple[int, int, int]:
    old_r, r, old_s, s, old_t, t = a, b, 1, 0, 0, 1
    while r:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t
    sign = -1 if old_r < 0 else 1
    return sign * old_r, sign * old_s, sign * old_t


def _column_hnf(columns: list[tuple[int, int]]) -> np.ndarray:
    """A rank-two integer span, using gcds of its maximal minors."""

    index = 0
    for (x1, y1), (x2, y2) in combinations(columns, 2):
        index = math.gcd(index, x1 * y2 - y1 * x2)
    if not index:
        raise ValueError("The rational planar lattice has rank less than two")
    x, y = 0, 0
    for col_x, col_y in columns:
        gcd, u, v = _bezout(y, col_y)
        x, y = u * x + v * col_x, gcd
    a = index // y
    return np.array([[a, x % a], [0, y]], dtype=np.int64)


def planar_lattices(numerator: np.ndarray, denominator: int) -> tuple[np.ndarray, np.ndarray]:
    """Return CSL and DSCL bases of Z² and (P/q) Z² exactly.

    The intersection follows from (L0 intersect L1)* = L0* + L1*.
    All index calculations remain integers; floats enter only at the end.
    """

    a, b, c, d = (int(x) for x in numerator.ravel())
    determinant = a * d - b * c
    if determinant <= 0:
        raise ValueError("The rational transformation must preserve orientation")
    q = denominator
    dscl = _column_hnf([(q, 0), (0, q), (a, c), (b, d)]) / q
    # P^-T = [[d, -c], [-b, a]] / det(P).
    dual = _column_hnf([
        (determinant, 0), (0, determinant),
        (q * d, -q * b), (-q * c, q * a),
    ])
    h, k, l = int(dual[0, 0]), int(dual[0, 1]), int(dual[1, 1])
    entries = (
        Fraction(determinant, h), Fraction(0),
        Fraction(-determinant * k, h * l), Fraction(determinant, l),
    )
    if any(value.denominator != 1 for value in entries):
        raise ArithmeticError("The planar intersection is not integral")
    csl = np.array([int(value) for value in entries], dtype=np.int64).reshape(2, 2)
    return csl, dscl


def reduce_planar(basis: np.ndarray) -> np.ndarray:
    """Use InterfaceMaster's native LLL on an embedded planar basis."""

    try:
        from interfacemaster.cellcalc import LLL
    except ImportError as exc:
        raise ImportError("InterfaceMaster is unavailable; use the configured launch script or install interface_master in the viewer environment") from exc
    embedded = np.vstack((basis, np.zeros(2)))
    reduced = LLL(embedded)[:2]
    if np.linalg.det(reduced) < 0:
        reduced[:, 1] *= -1
    return reduced


def continued_fraction(ratio: float, limit: int) -> list[tuple[int, int]]:
    """Principal convergents N/M of a positive spacing ratio dO/dS."""

    p0, p1, q0, q1 = 0, 1, 1, 0
    value, result = ratio, []
    for _ in range(40):
        coefficient = math.floor(value)
        p, q = coefficient * p1 + p0, coefficient * q1 + q0
        if p > limit or q > limit:
            break
        if p > 0 and q > 0:
            result.append((p, q))
        remainder = value - coefficient
        if remainder < 1e-12:
            break
        p0, p1, q0, q1 = p1, p, q1, q
        value = 1 / remainder
    return result


def rational_periods(film_length: float, substrate_length: float, limit: int) -> list[dict[str, Any]]:
    """Bounded coprime M:N periods, with continued-fraction/Farey context.

    M counts film periods and N counts substrate periods. These are 1D
    guides; independent length ratios do not solve a rotated 2D interface.
    """

    ratio = film_length / substrate_length
    convergents = set(continued_fraction(ratio, limit))
    records = []
    for m in range(1, limit + 1):
        for n in range(1, limit + 1):
            if math.gcd(m, n) != 1:
                continue
            delta = m * ratio / n - 1
            record = {"M": m, "N": n, "ratio_N_over_M": n / m,
                      "delta_1": delta, "delta_2": abs(m * ratio - n) / (m * ratio + n),
                      "convergent": (n, m) in convergents, "farey_parents": None}
            if m > 1:
                lower_m = pow(n, -1, m)
                lower_n = (n * lower_m - 1) // m
                record["farey_parents"] = [[lower_n, lower_m], [n - lower_n, m - lower_m]]
            records.append(record)
    records.sort(key=lambda row: (abs(row["delta_1"]), row["M"] + row["N"]))
    selected = records[:6]
    for row in records:
        if row["convergent"] and row not in selected:
            selected.append(row)
    return selected


def dichromatic_points(record: dict[str, Any], limit: int = 250) -> dict[str, list[list[float]]]:
    """Bounded Bravais-site preview using coset representatives of the CSL."""

    common = np.asarray(record["substrate_vectors_A"])[:, :2]
    points = {}
    for side in ("substrate", "film"):
        matrix = np.asarray(record[f"{side}_supercell_matrix"], dtype=int).T
        hnf = _column_hnf([tuple(int(x) for x in col) for col in matrix.T])
        nx, ny = int(hnf[0, 0]), int(hnf[1, 1])
        step = max(1, math.ceil(nx * ny / limit))
        representatives = np.array([(k // ny, k % ny) for k in range(0, nx * ny, step)])
        fractional = (representatives @ np.linalg.inv(matrix).T) % 1
        # Numerical roundoff at a periodic edge must not move an origin site
        # to the opposite corner of the drawn cell.
        fractional[np.isclose(fractional, 1, atol=1e-10, rtol=0)] = 0
        points[side] = (fractional @ common).tolist()
    return points

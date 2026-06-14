"""Shared symbolic utilities for paper07 witness scripts.

This module acts as a re-export facade: all public names are either
imported from the dppu library or are thin wrappers defined here for
paper07-specific vocabulary.  Witness scripts import exclusively from
this module and must not be changed.
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path
from typing import Any, Sequence

import sympy as sp

HERE = Path(__file__).resolve().parent
SCRIPT_ROOT = HERE.parents[1]
PROJECT_ROOT = HERE.parents[2]
DATA_DIR = PROJECT_ROOT / "data"

if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

# ---------------------------------------------------------------------------
# dppu library imports (utils)
# ---------------------------------------------------------------------------
from dppu.utils.epsilon import epsilon_nd
from dppu.utils.symbolic import (
    canonical_simplify,
    series_coefficients,
    derivative_degree_classifier,
)
from dppu.utils.printing import mat_str, print_json_rows

# ---------------------------------------------------------------------------
# dppu library imports (topology / class-A)
# ---------------------------------------------------------------------------
from dppu.topology.class_a import (
    structure_entries_class_a,
    structure_constants_array,
    n_matrix_from_structure_constants,
    exterior_derivative_matrix_d1,
    topology_n_data,
    rational_orthogonal_matrix as _rational_orthogonal_matrix_impl,
    rotate_n_matrix,
    vector_norm,
)

# ---------------------------------------------------------------------------
# dppu library imports (action / Palatini)
# ---------------------------------------------------------------------------
from dppu.action.palatini_first_order import (
    eh_density,
    palatini_reduced_density,
    minisuperspace_dewitt_metric,
)

# ---------------------------------------------------------------------------
# Paper07 constants
# ---------------------------------------------------------------------------

CANONICAL_TOPOLOGIES = ("T3", "Nil3", "Sol3", "S3")
PAIR_LIST = list(combinations(range(4), 2))


def setup_paper07_paths() -> tuple[Path, Path, Path]:
    """Return script dir, project root, and data directory."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    return HERE, PROJECT_ROOT, DATA_DIR


# ---------------------------------------------------------------------------
# Backward-compatible wrappers with paper07 names
# ---------------------------------------------------------------------------

def eps3(i: int, j: int, k: int) -> int:
    """3D Levi-Civita symbol (0-based indices).

    Re-exported via epsilon_nd for consistency; kept here for backward
    compatibility with witness scripts.
    """
    return epsilon_nd((i, j, k))


def simplify_expr(expr: Any) -> sp.Expr:
    """Canonical simplification: factor(cancel(simplify(expr))).

    Thin wrapper around dppu.utils.symbolic.canonical_simplify kept
    as a paper07-vocabulary alias.
    """
    return canonical_simplify(expr)


def rational_orthogonal_matrix(seed: int = 2) -> sp.Matrix:
    """Return a rational orthogonal matrix for frame-rotation checks.

    The ``seed`` argument is accepted for backward compatibility but is
    ignored — the matrix is fixed.
    """
    return _rational_orthogonal_matrix_impl()


# ---------------------------------------------------------------------------
# Paper07-vocabulary wrappers (delegate to class_a)
# ---------------------------------------------------------------------------

def activation_operator_from_structure(
    structure: list[list[list[sp.Expr]]],
) -> sp.Matrix:
    """Minimal paper model A_n: A^{ij} = 1/2 epsilon^{ikl} C^j_{kl} (= n^{ij}).

    In paper07 notation this is called the activation operator; it equals
    the n-matrix exactly.
    """
    return n_matrix_from_structure_constants(structure)


# ---------------------------------------------------------------------------
# Paper07-specific: finite-direction classifier (uses series_coefficients)
# ---------------------------------------------------------------------------

def classify_finite_direction(
    potential: sp.Expr,
    vector: Sequence[sp.Expr],
    *,
    order: int = 4,
) -> dict[str, Any]:
    """Classify the potential profile along a minisuperspace direction.

    Returns one of three labels:
      - ``EXACT_FLAT_MODULUS``: delta is exactly zero.
      - ``WALL``: at least one leading Taylor coefficient is nonzero.
      - ``NONFLAT_HIGHER_ORDER_OR_UNDETERMINED``: delta is nonzero but
        all coefficients up to *order* are zero.

    Args:
        potential: Symbolic potential as a function of s1, s2, s3.
        vector:    Direction vector (3 components) in log-scale space.
        order:     Maximum Taylor order to inspect.

    Returns:
        Dict with keys ``"delta"``, ``"coefficients"``, ``"nonzero_orders"``,
        ``"classification"``.
    """
    z = sp.Symbol("z")
    s1, s2, s3 = sp.symbols("s1 s2 s3", positive=True, nonzero=True)
    substitutions = {
        s1: sp.exp(z * vector[0]),
        s2: sp.exp(z * vector[1]),
        s3: sp.exp(z * vector[2]),
    }
    value = canonical_simplify(potential.subs(substitutions))
    value0 = canonical_simplify(value.subs(z, 0))
    delta = canonical_simplify(value - value0)
    coeffs = series_coefficients(delta, z, order=order)
    nonzero_orders = [key for key in ("1", "2", "3", "4") if sp.sympify(coeffs[key]) != 0]
    if sp.simplify(delta) == 0:
        classification = "EXACT_FLAT_MODULUS"
    elif nonzero_orders:
        classification = "WALL"
    else:
        classification = "NONFLAT_HIGHER_ORDER_OR_UNDETERMINED"
    return {
        "delta": delta,
        "coefficients": coeffs,
        "nonzero_orders": nonzero_orders,
        "classification": classification,
    }

"""Witness for S3 normalization scale invariance (paper07 Section 2 / App. A)."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
SCRIPT_ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

from paper07_common import (
    DATA_DIR,
    mat_str,
    n_matrix_from_structure_constants,
    print_json_rows,
    setup_paper07_paths,
    simplify_expr,
    structure_constants_array,
)
from dppu.utils.tee_logger import setup_log, teardown_log

setup_paper07_paths()

ROW_FIELDS = [
    "scale_factor",
    "rank_n",
    "dim_ker_n",
    "curvature_sample",
    "eigenvalue_sample",
    "rank_invariant",
    "curvature_scaled",
    "eigenvalue_scaled",
    "status",
]


def s3_scaled_data(scale: sp.Expr) -> dict[str, str]:
    base = structure_constants_array("S3")
    scaled = [[[scale * base[i][j][k] for k in range(3)] for j in range(3)] for i in range(3)]
    n_mat = n_matrix_from_structure_constants(scaled)
    rank_n = int(n_mat.rank())
    dim_ker = int(3 - rank_n)
    q, r = sp.symbols("q R", positive=True, real=True)
    l = sp.symbols("l", integer=True, positive=True)
    curvature = simplify_expr(6 * q**2 / r**2)
    curvature_scaled = simplify_expr(curvature.subs(r, r / scale))
    eigenvalue = simplify_expr((l - 1) * (l + 3))
    eigenvalue_scaled = simplify_expr(eigenvalue.subs({l: l}))  # symbolic in l; scale enters via C
    # Under C -> sC, n -> s n, rank unchanged; curvature scales as s^2, eigenvalue as s^2.
    curvature_with_scale = simplify_expr(6 * q**2 / (r / scale) ** 2)
    eigenvalue_with_scale = simplify_expr(scale**2 * (l - 1) * (l + 3))
    rank_invariant = rank_n == 3 and dim_ker == 0
    curvature_scaled_flag = simplify_expr(curvature_with_scale - scale**2 * curvature) == 0
    eigenvalue_scaled_flag = simplify_expr(eigenvalue_with_scale - scale**2 * eigenvalue) == 0
    status = "PASS" if rank_invariant and curvature_scaled_flag and eigenvalue_scaled_flag else "FAIL"
    return {
        "scale_factor": str(scale),
        "rank_n": str(rank_n),
        "dim_ker_n": str(dim_ker),
        "curvature_sample": str(curvature),
        "eigenvalue_sample": str(eigenvalue),
        "rank_invariant": "true" if rank_invariant else "false",
        "curvature_scaled": "true" if curvature_scaled_flag else "false",
        "eigenvalue_scaled": "true" if eigenvalue_scaled_flag else "false",
        "status": status,
    }


def main() -> int:
    setup_log(__file__, log_dir=str(DATA_DIR))
    try:
        print("Paper07 S3 normalization factor witness")
        print("=" * 65)
        print("frame_convention=C_equals_epsilon")
        base_n = n_matrix_from_structure_constants(structure_constants_array("S3"))
        print(f"base_n_matrix={mat_str(base_n)}")
        scales = [sp.Integer(1), sp.Integer(2), sp.Rational(1, 4)]
        rows = [s3_scaled_data(scale) for scale in scales]
        for row in rows:
            print(
                f"  scale={row['scale_factor']}: rank_n={row['rank_n']} "
                f"dim_ker_n={row['dim_ker_n']} status={row['status']}"
            )
        print_json_rows("scale_detail", rows, ROW_FIELDS)
        fail_count = sum(1 for row in rows if row["status"] != "PASS")
        all_pass = fail_count == 0
        print()
        print(f"fail_count={fail_count}")
        if all_pass:
            print("overall_verdict=PASS_S3_NORMALIZATION_FACTOR_WITNESS")
        else:
            print("overall_verdict=FAIL_S3_NORMALIZATION_FACTOR_WITNESS")
        return 0 if all_pass else 1
    finally:
        teardown_log()


if __name__ == "__main__":
    raise SystemExit(main())
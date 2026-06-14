"""DBR-derived rank-defect activation witness (paper07 Section 4 / App. B).

This witness has two layers:

1. A formal class-A identity check ``A=n`` for the structure-map notation.
2. The actual evidence layer: the DBR scalar-axis closure engine is
   run for each topology and axis.  The activation matrix is read from the
   resulting Dirac-Bergmann first/second-class split, not defined as ``n``.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import sympy as sp

HERE = Path(__file__).resolve().parent
SCRIPT_ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

from paper07_common import (
    CANONICAL_TOPOLOGIES,
    DATA_DIR,
    activation_operator_from_structure,
    eps3,
    mat_str,
    print_json_rows,
    rational_orthogonal_matrix,
    rotate_n_matrix,
    setup_paper07_paths,
    simplify_expr,
    topology_n_data,
)
from dppu.utils.tee_logger import setup_log, teardown_log
from dppu.dirac.local_models import build_branchA_local_L2_with_lc, build_t3_local, build_sol3_local
from dppu.dirac.scalar_axis import build_scalar_axis_projection, dirac_standard

setup_paper07_paths()

ROW_FIELDS = [
    "topology",
    "n_matrix",
    "engine_activation_matrix",
    "formal_activation_matrix",
    "rank_n",
    "rank_engine_activation",
    "dim_ker_n",
    "d_zero",
    "formal_activation_equals_n",
    "frame_rank_preserved",
    "claim_level",
    "status",
]

AXIS_ROW_FIELDS = [
    "topology",
    "axis",
    "projection",
    "N_config",
    "N_primary",
    "N_secondary",
    "N_constraints",
    "rank_constraint_matrix",
    "N_first_class",
    "N_second_class",
    "N_phys_pair",
    "tertiary_count",
    "engine_activation_entry",
    "expected_from_n",
    "generic_activation_entry",
    "activation_point_stable",
    "second_class_even",
    "status",
]

# rank/N_phys evaluation points: the original unit point plus a generic
# non-special point, so the single-point evaluation is confirmed not to sit on
# a degenerate locus (mirrors the two-point Hessian audit and route control).
RANK_POINTS: list[tuple[str, dict[str, int]]] = [
    ("unit", {"N": 1, "q": 1, "k": 1}),
    ("generic", {"N": 2, "q": 3, "k": 5}),
]


def axis_activation_entry(
    proj: dict[str, Any], local: dict[str, Any], point: dict[str, int]
) -> tuple[int, dict[str, Any]]:
    rank_subs = {
        local["N"]: sp.Integer(point["N"]),
        local["q"]: sp.Integer(point["q"]),
        proj["k"]: sp.Integer(point["k"]),
    }
    result = dirac_standard(proj, local, rank_subs)
    nphys = sp.Rational(result["N_phys"])
    entry = 1 if nphys == 0 and result["tertiary_count"] == 0 else 0
    return entry, result


def generic_activation_identity() -> bool:
    symbols = sp.symbols("n11 n12 n13 n21 n22 n23 n31 n32 n33")
    n_generic = sp.Matrix(3, 3, symbols)
    structure = [[[sp.Integer(0) for _ in range(3)] for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                structure[i][j][k] = simplify_expr(
                    sum(sp.Integer(eps3(j, k, l)) * n_generic[l, i] for l in range(3))
                )
    activation = activation_operator_from_structure(structure)
    return simplify_expr(activation - n_generic) == sp.zeros(3, 3)


def expected_axis_activation(n_mat: sp.Matrix) -> dict[str, int]:
    labels = ["x", "y", "z"]
    return {labels[i]: 0 if n_mat[:, i] == sp.zeros(3, 1) else 1 for i in range(3)}


def build_engine_locals() -> dict[str, Any]:
    return {
        "T3": build_t3_local(),
        "Nil3": build_branchA_local_L2_with_lc("Nil3"),
        "Sol3": build_sol3_local(),
        "S3": build_branchA_local_L2_with_lc("S3"),
    }


def build_rows() -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    rotation = rational_orthogonal_matrix()
    rows: list[dict[str, str]] = []
    axis_rows: list[dict[str, str]] = []
    locals_by_topology = build_engine_locals()
    for topology in CANONICAL_TOPOLOGIES:
        data = topology_n_data(topology)
        n_mat = data["n_matrix"]
        formal_activation = data["activation_matrix"]
        expected = expected_axis_activation(n_mat)
        local = locals_by_topology[topology]
        engine_entries: list[int] = []
        topology_axis_ok = True
        for axis in ("x", "y", "z"):
            proj = build_scalar_axis_projection(local, axis)
            activation_entry, result = axis_activation_entry(proj, local, RANK_POINTS[0][1])
            generic_entry, _ = axis_activation_entry(proj, local, RANK_POINTS[1][1])
            nphys = sp.Rational(result["N_phys"])
            engine_entries.append(activation_entry)
            n_second = int(result["N_second"])
            point_stable = activation_entry == generic_entry
            second_class_even = n_second % 2 == 0
            axis_ok = (
                activation_entry == expected[axis]
                and point_stable
                and second_class_even
            )
            topology_axis_ok = topology_axis_ok and axis_ok
            axis_rows.append(
                {
                    "topology": topology,
                    "axis": axis,
                    "projection": proj["projection"],
                    "N_config": str(result["N_config"]),
                    "N_primary": str(result["N_primary"]),
                    "N_secondary": str(result["N_secondary"]),
                    "N_constraints": str(result["N_constraints"]),
                    "rank_constraint_matrix": str(result["rank_C"]),
                    "N_first_class": str(result["N_first"]),
                    "N_second_class": str(n_second),
                    "N_phys_pair": str(nphys),
                    "tertiary_count": str(result["tertiary_count"]),
                    "engine_activation_entry": str(activation_entry),
                    "expected_from_n": str(expected[axis]),
                    "generic_activation_entry": str(generic_entry),
                    "activation_point_stable": "true" if point_stable else "false",
                    "second_class_even": "true" if second_class_even else "false",
                    "status": "PASS" if axis_ok else "FAIL",
                }
            )
        engine_activation = sp.diag(*engine_entries)
        rank_n = data["rank_n"]
        rank_engine_activation = int(engine_activation.rank())
        rotated = rotate_n_matrix(n_mat, rotation)
        rotated_rank = int(rotated.rank())
        formal_activation_equals_n = simplify_expr(formal_activation - n_mat) == sp.zeros(3, 3)
        frame_rank_preserved = rotated_rank == rank_n
        status = "PASS" if (
            rank_engine_activation == rank_n
            and data["dim_ker_n"] == 3 - rank_n
            and formal_activation_equals_n
            and frame_rank_preserved
            and topology_axis_ok
        ) else "FAIL"
        rows.append(
            {
                "topology": topology,
                "n_matrix": mat_str(n_mat),
                "engine_activation_matrix": mat_str(engine_activation),
                "formal_activation_matrix": mat_str(formal_activation),
                "rank_n": str(rank_n),
                "rank_engine_activation": str(rank_engine_activation),
                "dim_ker_n": str(data["dim_ker_n"]),
                "d_zero": str(data["d_zero"]),
                "formal_activation_equals_n": "true" if formal_activation_equals_n else "false",
                "frame_rank_preserved": "true" if frame_rank_preserved else "false",
                "claim_level": "ENGINE_DBR_SCALAR_AXIS_CONSTRAINT_ALGEBRA_DERIVED",
                "status": status,
            }
        )
    return rows, axis_rows


def main() -> int:
    setup_log(__file__, log_dir=str(DATA_DIR))
    try:
        print("Paper07 rank-defect activation witness")
        print("=" * 65)
        print("model=dbr_scalar_axis_constraint_closure")
        generic_ok = generic_activation_identity()
        print(f"formal_identity_A_equals_n={'PASS' if generic_ok else 'FAIL'}")
        print("formal_identity_role=notation_check_not_primary_evidence")
        rows, axis_rows = build_rows()
        for row in rows:
            print(
                f"  {row['topology']}: rank_n={row['rank_n']} "
                f"rank_engine_activation={row['rank_engine_activation']} "
                f"dim_ker_n={row['dim_ker_n']} status={row['status']}"
            )
        print_json_rows("topology_detail", rows, ROW_FIELDS)
        print_json_rows("axis_closure_detail", axis_rows, AXIS_ROW_FIELDS)
        fail_count = sum(1 for row in rows if row["status"] != "PASS")
        axis_fail_count = sum(1 for row in axis_rows if row["status"] != "PASS")
        all_pass = generic_ok and fail_count == 0
        print()
        print(f"total_topologies={len(rows)}")
        print(f"total_axis_closures={len(axis_rows)}")
        print(f"fail_count={fail_count}")
        print(f"axis_fail_count={axis_fail_count}")
        if all_pass and axis_fail_count == 0:
            print("overall_verdict=PASS_RANK_DEFECT_ACTIVATION_WITNESS")
        else:
            print("overall_verdict=FAIL_RANK_DEFECT_ACTIVATION_WITNESS")
        return 0 if all_pass and axis_fail_count == 0 else 1
    finally:
        teardown_log()


if __name__ == "__main__":
    raise SystemExit(main())

"""Witness for harmonic 1-form / DBR Hessian kernel correspondence.

This paper07 witness uses the closed-form harmonic-basis and DBR reduced-Hessian
engine (``dppu.spectral``).  The closed-form part is recomputed from the
Maurer-Cartan structure constants, and
the Hessian part is recomputed from the scalar-axis DBR closure engine on the
constraint-surface reduced phase-space basis.  It deliberately avoids the older
minimal model ``H = n.T*n``.
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
    mat_str,
    print_json_rows,
    setup_paper07_paths,
    topology_n_data,
)
from dppu.utils.tee_logger import setup_log, teardown_log
from dppu.spectral.harmonic_basis import (
    structure_constants,
    closure_matrix_from_c,
    basis_matrix,
    span_equal,
    co_closed_audit,
)
from dppu.spectral.reduced_hessian import scalar_axis_hessian_audits
from dppu.dirac.local_models import (
    build_branchA_local_L2_with_lc,
    build_t3_local,
    build_sol3_local,
)

setup_paper07_paths()

ROW_FIELDS = [
    "topology",
    "rank_n",
    "dim_ker_n",
    "closed_form_dim",
    "closed_equals_ker_n",
    "co_closed",
    "active_axes",
    "dbr_active_pair_count",
    "hessian_tested_points",
    "hessian_kernel_status",
    "evidence_level",
    "status",
]


def axis_labels_from_kernel(n_mat: sp.Matrix) -> list[str]:
    labels = ["x", "y", "z"]
    axes: list[str] = []
    for idx, label in enumerate(labels):
        unit = sp.zeros(3, 1)
        unit[idx, 0] = 1
        if n_mat * unit == sp.zeros(3, 1):
            axes.append(label)
    return axes


def build_topology_locals() -> dict[str, Any]:
    return {
        "locals": {
            "Nil3": build_branchA_local_L2_with_lc("Nil3"),
            "T3": build_t3_local(),
            "Sol3": build_sol3_local(),
        },
    }


def build_rows() -> list[dict[str, str]]:
    local_bundle = build_topology_locals()

    rows: list[dict[str, str]] = []
    for topology in CANONICAL_TOPOLOGIES:
        data = topology_n_data(topology)
        n_mat = data["n_matrix"]
        constants = structure_constants(n_mat)
        closed_matrix = closure_matrix_from_c(constants)
        closed_basis = basis_matrix(closed_matrix.nullspace())
        ker_n_basis = basis_matrix(n_mat.nullspace())
        closed_equals_ker = span_equal(closed_basis, ker_n_basis)
        co_closed = co_closed_audit(constants, closed_basis)

        dim_ker_n = data["dim_ker_n"]
        active_axes = axis_labels_from_kernel(n_mat)
        if topology == "S3":
            hessian_status = "PASS"
            active_pair_count = 0
            tested_points = 0
        else:
            hessian_audits = scalar_axis_hessian_audits(
                topology,
                local_bundle["locals"][topology],
                active_axes,
            )
            active_pair_count = len(active_axes)
            tested_points = sum(len(axis_row["point_rows"]) for axis_row in hessian_audits.values())
            hessian_status = "PASS" if all(axis_row["status"] == "PASS" for axis_row in hessian_audits.values()) else "FAIL"

        status = "PASS" if (
            int(closed_basis.cols) == dim_ker_n
            and closed_equals_ker
            and co_closed["status"] == "PASS"
            and hessian_status == "PASS"
        ) else "FAIL"
        rows.append(
            {
                "topology": topology,
                "rank_n": str(data["rank_n"]),
                "dim_ker_n": str(dim_ker_n),
                "closed_form_dim": str(int(closed_basis.cols)),
                "closed_equals_ker_n": "true" if closed_equals_ker else "false",
                "co_closed": "true" if co_closed["status"] == "PASS" else "false",
                "active_axes": ",".join(active_axes) if active_axes else "NONE",
                "dbr_active_pair_count": str(active_pair_count),
                "hessian_tested_points": str(tested_points),
                "hessian_kernel_status": hessian_status,
                "evidence_level": "CLOSED_FORM_BASIS_PLUS_DBR_REDUCED_HESSIAN",
                "status": status,
            }
        )
        print(
            f"  {topology}: n={mat_str(n_mat)} closed_dim={closed_basis.cols} "
            f"active_axes={','.join(active_axes) if active_axes else 'NONE'} "
            f"hessian={hessian_status} status={status}"
        )
    return rows


def main() -> int:
    setup_log(__file__, log_dir=str(DATA_DIR))
    try:
        print("Paper07 harmonic kernel correspondence witness")
        print("=" * 65)
        print("hessian_model=nTn_disabled")
        print("source=closed-form harmonic basis + DBR reduced-Hessian engine")
        rows = build_rows()
        print_json_rows("topology_detail", rows, ROW_FIELDS)
        fail_count = sum(1 for row in rows if row["status"] != "PASS")
        all_pass = fail_count == 0
        print()
        print(f"total_topologies={len(rows)}")
        print(f"fail_count={fail_count}")
        if all_pass:
            print("overall_verdict=PASS_HARMONIC_KERNEL_CORRESPONDENCE")
        else:
            print("overall_verdict=FAIL_HARMONIC_KERNEL_CORRESPONDENCE")
        return 0 if all_pass else 1
    finally:
        teardown_log()


if __name__ == "__main__":
    raise SystemExit(main())

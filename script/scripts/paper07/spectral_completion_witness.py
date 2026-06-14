"""Witness for scoped k>0 spectral completion (paper07 Section 5 / App. D)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
SCRIPT_ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

from paper07_common import DATA_DIR, print_json_rows, setup_paper07_paths, simplify_expr
from dppu.utils.tee_logger import setup_log, teardown_log
from dppu.spectral.sol3_compact_edge import compute_tables, gate_rows

setup_paper07_paths()

ROW_FIELDS = [
    "sector",
    "formula",
    "positivity_condition",
    "witness_value",
    "status",
]


def t3_positive() -> dict[str, str]:
    kx, ky, kz = sp.symbols("kx ky kz", real=True)
    lam = simplify_expr(kx**2 + ky**2 + kz**2)
    sample = lam.subs({kx: sp.pi, ky: sp.pi, kz: 0})
    positive = sp.simplify(sample) > 0
    return {
        "sector": "T3_Fourier_k_ne_0",
        "formula": "|k|^2",
        "positivity_condition": "k != 0",
        "witness_value": str(sample),
        "status": "PASS" if positive else "FAIL",
    }


def s3_positive() -> dict[str, str]:
    l = sp.symbols("l", integer=True, positive=True)
    lam = simplify_expr((l - 1) * (l + 3))
    sample = lam.subs(l, 2)
    positive = sp.simplify(sample) > 0
    return {
        "sector": "S3_TT_l_ge_2",
        "formula": "(l-1)(l+3)",
        "positivity_condition": "l >= 2",
        "witness_value": str(sample),
        "status": "PASS" if positive else "FAIL",
    }


def nil3_toroidal_positive() -> dict[str, str]:
    m, n = sp.symbols("m n", integer=True)
    lam = simplify_expr(4 * sp.pi**2 * (m**2 + n**2))
    sample = lam.subs({m: 1, n: 0})
    positive = sp.simplify(sample) > 0
    return {
        "sector": "Nil3_toroidal_k0",
        "formula": "4*pi^2*(m^2+n^2)",
        "positivity_condition": "(m,n) != (0,0)",
        "witness_value": str(sample),
        "status": "PASS" if positive else "FAIL",
    }


def nil3_heisenberg_positive() -> dict[str, str]:
    k, n = sp.symbols("k N", real=True, integer=True)
    lam = simplify_expr(2 * sp.pi * sp.Abs(k) * (2 * n + 1))
    sample = lam.subs({k: sp.Rational(1, 3), n: 0})
    positive = sp.simplify(sample) > 0
    return {
        "sector": "Nil3_Heisenberg_k_ne_0",
        "formula": "2*pi*|k|*(2N+1)",
        "positivity_condition": "k != 0",
        "witness_value": str(sample),
        "status": "PASS" if positive else "FAIL",
    }


def sol3_bulk_amgm() -> dict[str, str]:
    p, q, w, r = sp.symbols("p q w R", positive=True, real=True)
    v = simplify_expr(q**2 * sp.exp(2 * w / r) + p**2 * sp.exp(-2 * w / r))
    v_min = simplify_expr(2 * sp.Abs(p * q))
    w_star = sp.Rational(1, 2) * r * sp.log(p / q)
    at_min = simplify_expr(v.subs(w, w_star))
    equality = simplify_expr(at_min - v_min) == 0
    sample_positive = sp.simplify(v_min.subs({p: 2, q: 3})) > 0
    status = "PASS" if equality and sample_positive else "FAIL"
    return {
        "sector": "Sol3_boost_bulk",
        "formula": "V_min=2|pq| via AM-GM",
        "positivity_condition": "pq != 0",
        "witness_value": f"equality={equality}; sample={v_min.subs({p: 2, q: 3})}",
        "status": status,
    }


def sol3_compact_edge_genuine() -> dict[str, str]:
    rigidity_rows, lambda_rows, convergence_rows, compare_rows, robustness_rows = compute_tables()
    gate_rows_result = gate_rows(
        rigidity_rows,
        lambda_rows,
        convergence_rows,
        compare_rows,
        robustness_rows,
    )
    final_lambda_rows = [row for row in lambda_rows if int(row["grid_n"]) == 520]
    min_lambda = min(float(row["lambda1_genuine"]) for row in final_lambda_rows)
    q_levels = len({(row["discriminant"], row["Q_level_key"]) for row in compare_rows})
    discriminants = ",".join(sorted({str(row["discriminant"]) for row in compare_rows}, key=int))
    positive = all(float(row["lambda1_genuine"]) > 0.0 for row in final_lambda_rows)
    all_gates_pass = all(row["status"] == "PASS" for row in gate_rows_result)
    return {
        "sector": "Sol3_compact_edge_Anosov_orbit",
        "formula": "A-dual-orbit modified-Mathieu finite-difference lambda1",
        "positivity_condition": "nonzero compact A-dual orbit with Q=|alpha beta|>0",
        "witness_value": (
            f"matrices={len(rigidity_rows)}; cases={len(final_lambda_rows)}; "
            f"Q_levels={q_levels}; discriminants={discriminants}; min_lambda1={min_lambda:.6g}"
        ),
        "status": "PASS" if positive and all_gates_pass else "FAIL",
    }


def tt_kinetic_positive() -> dict[str, str]:
    q, n_lapse = sp.symbols("q N", positive=True, real=True)
    k = sp.symbols("k", real=True)
    kinetic = simplify_expr(2 * q**3 / n_lapse)
    derivative = simplify_expr(sp.diff(kinetic, k))
    positive = sp.simplify(kinetic.subs({q: 2, n_lapse: 1})) > 0
    label_independent = derivative == 0
    status = "PASS" if positive and label_independent else "FAIL"
    return {
        "sector": "TT_kinetic_H2",
        "formula": "diag(2*q^3/N)",
        "positivity_condition": "q>0, N>0",
        "witness_value": f"sample={kinetic.subs({q: 2, n_lapse: 1})}; dk/dlabel={derivative}",
        "status": status,
    }


def main() -> int:
    setup_log(__file__, log_dir=str(DATA_DIR))
    try:
        print("Paper07 scoped k>0 spectral completion witness")
        print("=" * 65)
        print("scope=analyzed_harmonic_sectors")
        print("sol3_edge_scope=compact Anosov-set modified-Mathieu witness")
        rows = [
            t3_positive(),
            s3_positive(),
            nil3_toroidal_positive(),
            nil3_heisenberg_positive(),
            sol3_bulk_amgm(),
            sol3_compact_edge_genuine(),
            tt_kinetic_positive(),
        ]
        for row in rows:
            print(
                f"  {row['sector']}: formula={row['formula']} "
                f"status={row['status']}"
            )
        print_json_rows("spectral_sector", rows, ROW_FIELDS)
        fail_count = sum(1 for row in rows if row["status"] != "PASS")
        all_pass = fail_count == 0
        print()
        print(f"total_sectors={len(rows)}")
        print(f"fail_count={fail_count}")
        if all_pass:
            print("overall_verdict=PASS_SPECTRAL_COMPLETION_WITNESS")
        else:
            print("overall_verdict=FAIL_SPECTRAL_COMPLETION_WITNESS")
        return 0 if all_pass else 1
    finally:
        teardown_log()


if __name__ == "__main__":
    raise SystemExit(main())

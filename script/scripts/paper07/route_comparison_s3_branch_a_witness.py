"""Route-comparison witness from direct recomputation of the bundled DBR/RBD engine.

No intermediate CSV is used as evidence.  This script uses the bundled DBR
constraint-closure and RBD route-variant engine (``dppu.dirac``) and recomputes
the compared route variants from the same Branch-A local first-order Lagrangian:

* DBR_KEEP_ALL: keep the independent connection through Dirac closure.
* RBD_SPATIAL_OMEGA_ROUTH: Routh-eliminate spatial omega first, then Dirac.
* RBD_FULL_K_ROUTH: Routh-eliminate all K first, then Dirac.
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

from paper07_common import DATA_DIR, print_json_rows, setup_paper07_paths
from dppu.utils.tee_logger import setup_log, teardown_log
from dppu.dirac.constraint_closure import (
    build_branchA_local_L2,
    project_lagrangian,
    build_dirac_for_projection,
)
from dppu.dirac.scalar_axis import build_scalar_axis_projection, dirac_standard
from dppu.dirac.route_variants import routh_eliminate, project_reduced_dirac

setup_paper07_paths()

L_SECTORS = ("l0", "l1", "lge2")

ROUTE_FIELDS = [
    "route_id",
    "definition",
    "l0",
    "l1",
    "lge2",
    "source_stage",
    "status",
]

DETAIL_FIELDS = [
    "route_id",
    "l_sector",
    "sample_l",
    "N_config",
    "N_constraints",
    "rank_constraint_matrix",
    "N_first_class",
    "N_second_class",
    "N_phys_pair",
    "calculation_path",
    "status",
]

COMM_FIELDS = [
    "route_A",
    "route_B",
    "Nphys_A_by_l",
    "Nphys_B_by_l",
    "commutator",
    "notes",
    "status",
]

SIGN_FIELDS = [
    "route_id",
    "l_sector",
    "surviving_pair",
    "H_PP",
    "sign",
    "calculation",
    "status",
]

CONTROL_FIELDS = [
    "control_id",
    "axis",
    "rank_point",
    "N_config",
    "N_constraints",
    "rank_constraint_matrix",
    "N_phys_pair",
    "calculation_path",
    "status",
]


def lkey(label: str) -> str:
    return {"l=0": "l0", "l=1": "l1", "l>=2": "lge2"}[label]


def compute_route_data() -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    local = build_branchA_local_L2()
    routh_spatial = routh_eliminate(local, also_K0=False)
    routh_full_k = routh_eliminate(local, also_K0=True)
    if routh_spatial.get("status") != "REDUCED":
        raise RuntimeError(f"spatial omega Routh failed: {routh_spatial}")
    if routh_full_k.get("status") != "REDUCED":
        raise RuntimeError(f"full K Routh failed: {routh_full_k}")

    rows: list[dict[str, str]] = []
    sign_rows: list[dict[str, str]] = []
    for l_sector, sample_l in (("l=0", 0), ("l=1", 1), ("l>=2", 2)):
        proj = project_lagrangian(local, sample_l, l_sector)
        dbr = build_dirac_for_projection(proj, local)
        spatial = project_reduced_dirac(local, routh_spatial, sample_l)
        full_k = project_reduced_dirac(local, routh_full_k, sample_l)
        if sample_l == 0:
            sign_rows.extend(
                [
                    l0_hpp_sign_row("RBD_SPATIAL_OMEGA_ROUTH", spatial),
                    l0_hpp_sign_row("RBD_FULL_K_ROUTH", full_k),
                ]
            )
        rows.extend(
            [
                detail_row(
                    "DBR_KEEP_ALL",
                    l_sector,
                    sample_l,
                    dbr,
                    "direct full DBR closure; no connection Routh elimination",
                ),
                detail_row(
                    "RBD_SPATIAL_OMEGA_ROUTH",
                    l_sector,
                    sample_l,
                    spatial,
                    "direct spatial-omega Routh, then scalar projection and Dirac",
                ),
                detail_row(
                    "RBD_FULL_K_ROUTH",
                    l_sector,
                    sample_l,
                    full_k,
                    "direct full-K Routh, then scalar projection and Dirac",
                ),
            ]
        )
    control_rows = scalar_axis_control_rows(local)
    return rows, sign_rows, control_rows


def scalar_axis_control_rows(local: dict[str, Any]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    rank_points = [
        ("unit", {local["N"]: sp.S.One, local["q"]: sp.S.One}),
        ("generic", {local["N"]: sp.Integer(2), local["q"]: sp.Integer(3)}),
    ]
    for axis in ("x", "y", "z"):
        proj = build_scalar_axis_projection(local, axis)
        for point_id, base_subs in rank_points:
            rank_subs = dict(base_subs)
            rank_subs[proj["k"]] = sp.Integer(1) if point_id == "unit" else sp.Integer(5)
            result = dirac_standard(proj, local, rank_subs)
            rows.append(
                {
                    "control_id": "SCALAR_AXIS_EXPAND_POISSON_DBR_CONTROL",
                    "axis": axis,
                    "rank_point": point_id,
                    "N_config": str(result.get("N_config", "")),
                    "N_constraints": str(result.get("N_constraints", "")),
                    "rank_constraint_matrix": str(result.get("rank_C", "")),
                    "N_phys_pair": str(result.get("N_phys", "")),
                    "calculation_path": "unreduced Branch-A local L through dppu.dirac.scalar_axis.dirac_standard; expand-based Poisson; two rank points",
                    "status": "PASS" if str(result.get("N_phys", "")) == "0" else "FAIL",
                }
            )
    return rows


def detail_row(route_id: str, l_sector: str, sample_l: int, result: dict[str, Any], path: str) -> dict[str, str]:
    return {
        "route_id": route_id,
        "l_sector": l_sector,
        "sample_l": str(sample_l),
        "N_config": str(result.get("N_config", "")),
        "N_constraints": str(result.get("N_constraints", "")),
        "rank_constraint_matrix": str(result.get("rank_C", "")),
        "N_first_class": str(result.get("N_first", "")),
        "N_second_class": str(result.get("N_second", "")),
        "N_phys_pair": str(result.get("N_phys", "")),
        "calculation_path": path,
        "status": "COMPUTED" if result.get("N_phys") is not None else str(result.get("status", "FAIL")),
    }


def l0_hpp_sign_row(route_id: str, result: dict[str, Any]) -> dict[str, str]:
    momentum = sp.Symbol("P_s_l0")
    lapse_pert = sp.Symbol("n_l0")
    h_c = sp.sympify(result["H_c"])
    h_pp = sp.simplify(sp.diff(h_c, momentum, momentum).subs({lapse_pert: 0}))
    if h_pp < 0:
        sign = "NEGATIVE"
    elif h_pp > 0:
        sign = "POSITIVE"
    elif h_pp == 0:
        sign = "ZERO"
    else:
        sign = "UNDETERMINED"
    return {
        "route_id": route_id,
        "l_sector": "l=0",
        "surviving_pair": "(s_l0, P_s_l0)",
        "H_PP": str(h_pp),
        "sign": sign,
        "calculation": "direct second derivative d^2 H_c / d P_s_l0^2 at n_l0=0 from recomputed RBD Hamiltonian",
        "status": "PASS" if sign == "NEGATIVE" and result.get("N_phys") == 1 else "FAIL",
    }


def route_summary_rows(detail_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_route: dict[str, dict[str, str]] = {}
    definitions = {
        "DBR_KEEP_ALL": "Dirac-before-Routh: keep independent connection through DBR closure",
        "RBD_SPATIAL_OMEGA_ROUTH": "Routh-before-Dirac: eliminate spatial omega first",
        "RBD_FULL_K_ROUTH": "Routh-before-Dirac: eliminate all K first",
    }
    for row in detail_rows:
        route = by_route.setdefault(
            row["route_id"],
            {
                "route_id": row["route_id"],
                "definition": definitions[row["route_id"]],
                "l0": "MISSING",
                "l1": "MISSING",
                "lge2": "MISSING",
                "source_stage": "direct DBR/RBD engine recomputation",
                "status": "PASS",
            },
        )
        route[lkey(row["l_sector"])] = row["N_phys_pair"]
        if row["status"] != "COMPUTED":
            route["status"] = "FAIL"
    return [by_route[key] for key in ("DBR_KEEP_ALL", "RBD_SPATIAL_OMEGA_ROUTH", "RBD_FULL_K_ROUTH")]


def commutator(a: dict[str, str], b: dict[str, str]) -> tuple[str, str]:
    va = [a[key] for key in L_SECTORS]
    vb = [b[key] for key in L_SECTORS]
    if va != vb:
        diff = ", ".join(f"{key}:{a[key]}vs{b[key]}" for key in L_SECTORS if a[key] != b[key])
        return "NONCOMMUTES", f"Nphys differ ({diff})"
    return "COMMUTES", "Nphys equal at all l sectors"


def commutator_rows(routes: list[dict[str, str]]) -> list[dict[str, str]]:
    by_id = {row["route_id"]: row for row in routes}
    pairs = [
        ("DBR_KEEP_ALL", "RBD_SPATIAL_OMEGA_ROUTH"),
        ("DBR_KEEP_ALL", "RBD_FULL_K_ROUTH"),
        ("RBD_SPATIAL_OMEGA_ROUTH", "RBD_FULL_K_ROUTH"),
    ]
    rows: list[dict[str, str]] = []
    for left, right in pairs:
        verdict, note = commutator(by_id[left], by_id[right])
        rows.append(
            {
                "route_A": left,
                "route_B": right,
                "Nphys_A_by_l": "/".join(by_id[left][key] for key in L_SECTORS),
                "Nphys_B_by_l": "/".join(by_id[right][key] for key in L_SECTORS),
                "commutator": verdict,
                "notes": note,
                "status": "PASS",
            }
        )
    return rows


def main() -> int:
    setup_log(__file__, log_dir=str(DATA_DIR))
    try:
        print("Paper07 route comparison S3 Branch-A witness")
        print("=" * 65)
        print("scope=direct DBR/RBD engine recomputation")
        print("physical_instability_claim=false")
        print("toy_model_removed=true")
        print("intermediate_csv_evidence=false")
        details, sign_rows, control_rows = compute_route_data()
        routes = route_summary_rows(details)
        comms = commutator_rows(routes)
        for row in comms:
            print(f"  Comm({row['route_A']},{row['route_B']})={row['commutator']} [{row['notes']}]")
        print_json_rows("route_detail", details, DETAIL_FIELDS)
        print_json_rows("route_l0_sign", sign_rows, SIGN_FIELDS)
        print_json_rows("route_control", control_rows, CONTROL_FIELDS)
        print_json_rows("route_nphys", routes, ROUTE_FIELDS)
        print_json_rows("route_commutator", comms, COMM_FIELDS)
        checks = {
            "DBR_recomputed_all_zero": next(row for row in routes if row["route_id"] == "DBR_KEEP_ALL")["l0"] == "0"
            and next(row for row in routes if row["route_id"] == "DBR_KEEP_ALL")["l1"] == "0"
            and next(row for row in routes if row["route_id"] == "DBR_KEEP_ALL")["lge2"] == "0",
            "RBD_spatial_noncommutes_with_DBR": any(
                row["route_A"] == "DBR_KEEP_ALL"
                and row["route_B"] == "RBD_SPATIAL_OMEGA_ROUTH"
                and row["commutator"] == "NONCOMMUTES"
                for row in comms
            ),
            "RBD_fullK_noncommutes_with_DBR": any(
                row["route_A"] == "DBR_KEEP_ALL"
                and row["route_B"] == "RBD_FULL_K_ROUTH"
                and row["commutator"] == "NONCOMMUTES"
                for row in comms
            ),
            "RBD_l0_HPP_negative": all(row["status"] == "PASS" for row in sign_rows),
            "scalar_axis_independent_control_all_zero": all(row["status"] == "PASS" for row in control_rows),
        }
        for key, value in checks.items():
            print(f"{key}={value}")
        all_pass = (
            all(checks.values())
            and all(row["status"] == "PASS" for row in routes)
            and all(row["status"] == "PASS" for row in control_rows)
        )
        print()
        print(f"fail_count={0 if all_pass else 1}")
        if all_pass:
            print("overall_verdict=PASS_ROUTE_COMPARISON_S3_BRANCH_A_WITNESS")
        else:
            print("overall_verdict=FAIL_ROUTE_COMPARISON_S3_BRANCH_A_WITNESS")
        return 0 if all_pass else 1
    finally:
        teardown_log()


if __name__ == "__main__":
    raise SystemExit(main())

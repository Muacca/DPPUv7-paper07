"""Witness for torsion auxiliary shell and field promotion (paper07 Section 8 / App. G)."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Callable, Sequence

import sympy as sp
from sympy import Matrix, S, cancel
from sympy.tensor.array import MutableDenseNDimArray

HERE = Path(__file__).resolve().parent
SCRIPT_ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

from paper07_common import (
    CANONICAL_TOPOLOGIES,
    DATA_DIR,
    print_json_rows,
    setup_paper07_paths,
)
from dppu.action.reduced_sector import (
    MODES,
    derive_static_branch_function,
    solve_auxiliary_shell,
)
from dppu.curvature.ricci import compute_ricci_scalar
from dppu.curvature.riemann import compute_riemann_tensor
from dppu.engine.contortion import compute_contortion
from dppu.engine.ec_connection import compute_ec_connection
from dppu.engine.levi_civita import compute_christoffel_frame
from dppu.forms.exterior_algebra import (
    Form,
    add_forms,
    basis,
    clean_form,
    coefficient,
    scale_form,
    wedge,
)
from dppu.torsion.ansatz import construct_torsion_tensor
from dppu.torsion.mode import Mode
from dppu.utils.epsilon import epsilon_3d
from dppu.utils.tee_logger import setup_log, teardown_log

setup_paper07_paths()

DIM = 4
METRIC = Matrix.diag(-1, 1, 1, 1)
METRIC_INV = Matrix.diag(-1, 1, 1, 1)

A1S, A2S, A3S = sp.symbols("a1 a2 a3", positive=True)
NS = sp.symbols("N", positive=True)
ETAS = sp.symbols("eta", real=True)
VS = sp.symbols("V", real=True)
A1D, A2D, A3D = sp.symbols("a1d a2d a3d", real=True)
ND = sp.symbols("Nd", real=True)
ETAD = sp.symbols("etad", real=True)
VD = sp.symbols("Vd", real=True)
DETA = list(sp.symbols("deta1 deta2 deta3", real=True))
DV = list(sp.symbols("dV1 dV2 dV3", real=True))
DERIV_SYMS_TIME = [ETAD, VD]
DERIV_SYMS_SPACE = DETA + DV

ROW_FIELDS = [
    "topology",
    "mode",
    "max_degree_eta_V",
    "shell_eta_zero",
    "shell_V_zero",
    "quadratic_derivative_absent",
    "native_derivative_rows",
    "native_gate_source",
    "status",
]

NATIVE_FIELDS = [
    "gate_id",
    "check",
    "computed",
    "expected",
    "status",
]

DERIVATIVE_FIELDS = [
    "topology",
    "mode",
    "sector",
    "poly_residual_zero",
    "max_deriv_degree",
    "classification",
    "linear_coefficients",
    "quad_nonzero_count",
    "quad_nonzero_monomials",
    "time_quad_present",
    "space_quad_present",
    "mixed_ts_quad_present",
    "higher_order_nonzero",
    "ibp_residual_zero",
    "note",
    "status",
]


def max_total_degree(expr: sp.Expr, eta: sp.Symbol, v_var: sp.Symbol) -> int:
    lam = sp.Symbol("lam")
    expanded = sp.expand(expr.subs({eta: lam * eta, v_var: lam * v_var}))
    if expanded == 0:
        return 0
    poly = sp.Poly(expanded, lam)
    return max(sum(degree) for degree, _ in poly.terms())


def compact(expr: Any, max_len: int = 140) -> str:
    text = str(cancel(expr))
    return text if len(text) <= max_len else text[: max_len - 3] + "..."


def mode_object(mode: str) -> Mode:
    return Mode[mode] if isinstance(mode, str) else mode


def build_structure_constants(
    topology: str,
    a1: Any,
    a2: Any,
    a3: Any,
    lapse: Any,
    a1d: Any,
    a2d: Any,
    a3d: Any,
) -> MutableDenseNDimArray:
    c_arr = MutableDenseNDimArray.zeros(DIM, DIM, DIM)
    hubble = [a1d / (lapse * a1), a2d / (lapse * a2), a3d / (lapse * a3)]
    for i in range(1, 4):
        c_arr[i, 0, i] = hubble[i - 1]
        c_arr[i, i, 0] = -hubble[i - 1]

    scale = [None, a1, a2, a3]
    n_by_topology = {
        "T3": [0, 0, 0],
        "Nil3": [0, 0, 1],
        "Sol3": [0, 1, -1],
        "S3": [1, 1, 1],
    }
    n_diag = n_by_topology[topology]
    for i in range(1, 4):
        for j in range(1, 4):
            for k in range(1, 4):
                eps_val = epsilon_3d(i - 1, j - 1, k - 1)
                if eps_val != 0 and n_diag[k - 1] != 0:
                    c_arr[k, i, j] = cancel(eps_val * n_diag[k - 1] * scale[k] / (scale[i] * scale[j]))
    return c_arr


def build_torsion_anisotropic(mode: str, a1: Any, a2: Any, a3: Any) -> MutableDenseNDimArray:
    r_eff = cancel((a1 * a2 * a3) ** sp.Rational(1, 3))
    return construct_torsion_tensor(
        mode_object(mode),
        r_eff,
        ETAS,
        VS,
        METRIC,
        dim=DIM,
        frame_convention="lz_native",
        signature="lorentzian",
    )


def field_frame_deriv() -> Callable[[int, Any], Any]:
    a1dd, a2dd, a3dd = sp.symbols("a1dd a2dd a3dd", real=True)
    ndd = sp.symbols("Ndd", real=True)
    etadd, vdd = sp.symbols("etadd Vdd", real=True)
    time_pairs = [
        (A1S, A1D),
        (A2S, A2D),
        (A3S, A3D),
        (NS, ND),
        (ETAS, ETAD),
        (VS, VD),
        (A1D, a1dd),
        (A2D, a2dd),
        (A3D, a3dd),
        (ND, ndd),
        (ETAD, etadd),
        (VD, vdd),
    ]

    def _deriv(mu: int, expr: Any) -> Any:
        if expr == S.Zero:
            return S.Zero
        if mu == 0:
            result = S.Zero
            for sym, dsym in time_pairs:
                dpart = sp.diff(expr, sym)
                if dpart != S.Zero:
                    result += dpart * dsym
            return cancel(result / NS)
        result = S.Zero
        eta_part = sp.diff(expr, ETAS)
        if eta_part != S.Zero:
            result += eta_part * DETA[mu - 1]
        v_part = sp.diff(expr, VS)
        if v_part != S.Zero:
            result += v_part * DV[mu - 1]
        return cancel(result)

    return _deriv


def de_basis_forms(c_arr: MutableDenseNDimArray) -> tuple[Form, ...]:
    forms = []
    for a in range(DIM):
        item: Form = {}
        for b in range(DIM):
            for c in range(DIM):
                val = c_arr[a, b, c] / 2
                if val != S.Zero:
                    item = add_forms(item, scale_form(val, wedge(basis(b), basis(c))))
        forms.append(item)
    return tuple(forms)


def exterior_derivative_field(
    form: Form,
    c_arr: MutableDenseNDimArray,
    frame_deriv: Callable[[int, Any], Any],
) -> Form:
    de = de_basis_forms(c_arr)
    out: Form = {}
    for idx, coeff in form.items():
        if coeff != S.Zero:
            for a in range(DIM):
                df_a = frame_deriv(a, coeff)
                if df_a != S.Zero:
                    out = add_forms(out, wedge({(a,): df_a}, {idx: S.One}))
        for pos, a in enumerate(idx):
            prefix: Form = {idx[:pos]: S.One}
            suffix: Form = {idx[pos + 1:]: S.One}
            term = wedge(wedge(prefix, de[a]), suffix)
            if pos % 2:
                term = scale_form(-S.One, term)
            out = add_forms(out, scale_form(coeff, term))
    return clean_form(out)


def d_scalar_form(expr: Any, frame_deriv: Callable[[int, Any], Any]) -> Form:
    out: Form = {}
    for a in range(DIM):
        value = frame_deriv(a, expr)
        if value != S.Zero:
            out[(a,)] = value
    return clean_form(out)


def q_nieh_yan_form(tensor_low: MutableDenseNDimArray) -> Form:
    qform: Form = {}
    for a in range(DIM):
        torsion_form: Form = {}
        for b in range(DIM):
            for c in range(DIM):
                val = tensor_low[a, b, c] / 2
                if val != S.Zero:
                    torsion_form = add_forms(torsion_form, scale_form(val, wedge(basis(b), basis(c))))
        qform = add_forms(qform, wedge(basis(a), torsion_form))
    return clean_form(qform)


def classify_derivative_structure(expr: Any, max_probe_order: int = 4) -> dict[str, Any]:
    symbols = DERIV_SYMS_TIME + DERIV_SYMS_SPACE
    lam = sp.symbols("lam_deriv", real=True)
    expr_lam = expr.subs({sym: lam * sym for sym in symbols}, simultaneous=True)
    parts: dict[int, Any] = {}
    work = expr_lam
    for order in range(max_probe_order + 1):
        parts[order] = cancel(work.subs(lam, 0))
        work = sp.diff(work, lam) / (order + 1)
    residual = cancel(sp.together(expr_lam - sum(parts[i] * lam**i for i in range(max_probe_order + 1))))
    poly_ok = residual == S.Zero
    max_degree = max((order for order, part in parts.items() if part != S.Zero), default=0)
    linear = {
        str(sym): cancel(sp.diff(parts.get(1, S.Zero), sym))
        for sym in symbols
        if cancel(sp.diff(parts.get(1, S.Zero), sym)) != S.Zero
    }
    quad: dict[str, Any] = {}
    quad_part = parts.get(2, S.Zero)
    for i, left in enumerate(symbols):
        for right in symbols[i:]:
            coeff = cancel(sp.diff(quad_part, left, 2) / 2) if left == right else cancel(sp.diff(sp.diff(quad_part, left), right))
            if coeff != S.Zero:
                quad[f"{left}*{right}"] = coeff
    higher = {order: parts[order] for order in range(3, max_probe_order + 1) if parts[order] != S.Zero}
    time_names = {str(sym) for sym in DERIV_SYMS_TIME}
    space_names = {str(sym) for sym in DERIV_SYMS_SPACE}

    def names(key: str) -> tuple[str, str]:
        left, right = key.split("*")
        return left, right

    time_quad = any(names(key)[0] in time_names and names(key)[1] in time_names for key in quad)
    space_quad = any(names(key)[0] in space_names and names(key)[1] in space_names for key in quad)
    mixed_ts = any((names(key)[0] in time_names) != (names(key)[1] in time_names) for key in quad)
    if max_degree == 0:
        classification = "absent"
    elif max_degree == 1:
        classification = "linear_current"
    else:
        classification = "quadratic_kinetic"
    return {
        "poly_residual_zero": poly_ok,
        "max_degree": max_degree,
        "linear": linear,
        "quad": quad,
        "higher": higher,
        "time_quad_present": time_quad,
        "space_quad_present": space_quad,
        "mixed_ts_quad_present": mixed_ts,
        "classification": classification,
    }


def compute_delta_r_field(topology: str, mode: str) -> tuple[Any, str]:
    c_arr = build_structure_constants(topology, A1S, A2S, A3S, NS, A1D, A2D, A3D)
    frame_deriv = field_frame_deriv()
    tensor_low = build_torsion_anisotropic(mode, A1S, A2S, A3S)
    contortion = compute_contortion(tensor_low, DIM, metric=METRIC)
    gamma_lc = compute_christoffel_frame(c_arr, DIM, metric=METRIC, metric_inv=METRIC_INV)
    gamma_ec = compute_ec_connection(gamma_lc, contortion, DIM)
    r_lc = compute_riemann_tensor(gamma_lc, c_arr, DIM, frame_deriv=frame_deriv)
    r_ec = compute_riemann_tensor(gamma_ec, c_arr, DIM, frame_deriv=frame_deriv)
    delta_r = cancel(compute_ricci_scalar(r_ec, DIM, metric_inv=METRIC_INV) - compute_ricci_scalar(r_lc, DIM, metric_inv=METRIC_INV))
    nonzero_k = sum(1 for a in range(DIM) for b in range(DIM) for c in range(DIM) if contortion[a, b, c] != S.Zero)
    return delta_r, f"DeltaR via bundled Riemann engine; contortion_nonzero={nonzero_k}"


def build_ny_sector_forms(topology: str, mode: str) -> dict[str, Any]:
    c_arr = build_structure_constants(topology, A1S, A2S, A3S, NS, A1D, A2D, A3D)
    frame_deriv = field_frame_deriv()
    tensor_low = build_torsion_anisotropic(mode, A1S, A2S, A3S)
    qform = q_nieh_yan_form(tensor_low)
    d_q = exterior_derivative_field(qform, c_arr, frame_deriv)
    d_q_coeff = cancel(coefficient(d_q, (0, 1, 2, 3)))
    d_eta = d_scalar_form(ETAS, frame_deriv)
    d_eta_q_coeff = cancel(coefficient(wedge(d_eta, qform), (0, 1, 2, 3)))
    d_eta_q = exterior_derivative_field(scale_form(ETAS, qform), c_arr, frame_deriv)
    ibp_residual = cancel(coefficient(d_eta_q, (0, 1, 2, 3)) - (d_eta_q_coeff + ETAS * d_q_coeff))
    return {
        "dQ_coeff": d_q_coeff,
        "dEtaQ_coeff": d_eta_q_coeff,
        "eta_dQ_coeff": cancel(ETAS * d_q_coeff),
        "ibp_residual": ibp_residual,
    }


def background_subs(topology: str) -> dict[Any, Any]:
    q = sp.symbols("q", positive=True)
    return {
        "T3": {A1S: sp.S.One, A2S: sp.S.One, A3S: sp.S.One},
        "Nil3": {A1S: sp.S.One, A2S: sp.S.One, A3S: sp.S.One},
        "Sol3": {A1S: sp.S.One, A2S: sp.S.One, A3S: sp.S.One},
        "S3": {A1S: q, A2S: q, A3S: q},
    }[topology]


def derivative_row(topology: str, mode: str, sector: str, expr: Any, note: str, ibp_zero: bool) -> dict[str, str]:
    volume = A1S * A2S * A3S
    lagrangian = cancel(volume * NS * expr).subs(background_subs(topology))
    data = classify_derivative_structure(cancel(lagrangian))
    quad = data["quad"]
    status = "PASS" if (
        data["poly_residual_zero"]
        and not data["time_quad_present"]
        and not data["space_quad_present"]
        and not data["mixed_ts_quad_present"]
        and not data["higher"]
        and len(quad) == 0
    ) else "FAIL"
    return {
        "topology": topology,
        "mode": mode,
        "sector": sector,
        "poly_residual_zero": "yes" if data["poly_residual_zero"] else "no",
        "max_deriv_degree": str(data["max_degree"]),
        "classification": data["classification"],
        "linear_coefficients": "; ".join(f"{key}:{compact(value, 60)}" for key, value in sorted(data["linear"].items())) or "none",
        "quad_nonzero_count": str(len(quad)),
        "quad_nonzero_monomials": "; ".join(f"{key}:{compact(value, 60)}" for key, value in sorted(quad.items())) or "none",
        "time_quad_present": "yes" if data["time_quad_present"] else "no",
        "space_quad_present": "yes" if data["space_quad_present"] else "no",
        "mixed_ts_quad_present": "yes" if data["mixed_ts_quad_present"] else "no",
        "higher_order_nonzero": "yes" if data["higher"] else "no",
        "ibp_residual_zero": "yes" if ibp_zero else "no",
        "note": note,
        "status": status,
    }


def native_derivative_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for topology in CANONICAL_TOPOLOGIES:
        for mode in ("AX", "VT", "MX"):
            delta_r, delta_note = compute_delta_r_field(topology, mode)
            ny = build_ny_sector_forms(topology, mode)
            rows.append(derivative_row(topology, mode, "EC", delta_r, delta_note, ny["ibp_residual"] == S.Zero))
            rows.append(derivative_row(topology, mode, "NY_dQNY", ny["dQ_coeff"], "NY-MAIN dQ_NY coefficient from bundled form engine", ny["ibp_residual"] == S.Zero))
            rows.append(derivative_row(topology, mode, "NY_dEta_wedge_QNY", ny["dEtaQ_coeff"], "IBP current candidate d eta wedge Q_NY", ny["ibp_residual"] == S.Zero))
            rows.append(derivative_row(topology, mode, "NY_eta_dQNY", ny["eta_dQ_coeff"], "eta*dQ_NY coupling with in-process IBP residual", ny["ibp_residual"] == S.Zero))
    return rows


def detector_control_rows(derivative_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    phi, phidot, grad = sp.symbols("phi phidot grad", real=True)
    free_scalar = sp.Rational(1, 2) * phidot**2 - sp.Rational(1, 2) * grad**2
    free_data = classify_derivative_structure(free_scalar.subs({phidot: ETAD, grad: DETA[0]}))
    linear_current = classify_derivative_structure(ETAD + 2 * VD)
    nc_full_pde = "NOT_TRIGGERED"
    all_native_rows_pass = all(row["status"] == "PASS" for row in derivative_rows)
    any_linear = any(row["classification"] == "linear_current" for row in derivative_rows)
    all_spatial_linear_zero = all(
        all(f"deta{i}:" not in row["linear_coefficients"] for i in range(1, 4))
        and all(f"dV{i}:" not in row["linear_coefficients"] for i in range(1, 4))
        for row in derivative_rows
    )
    return [
        {
            "gate_id": "NG-01",
            "check": "in-process native derivative rows",
            "computed": f"rows={len(derivative_rows)}; pass={sum(1 for row in derivative_rows if row['status'] == 'PASS')}",
            "expected": "all rows PASS; no quadratic derivative monomials",
            "status": "PASS" if all_native_rows_pass else "FAIL",
        },
        {
            "gate_id": "PC-01",
            "check": "positive detector fires on free scalar kinetic/gradient",
            "computed": f"class={free_data['classification']}; quad={len(free_data['quad'])}",
            "expected": "quadratic_kinetic with nonzero quadratic monomials",
            "status": "PASS" if free_data["classification"] == "quadratic_kinetic" and len(free_data["quad"]) == 2 else "FAIL",
        },
        {
            "gate_id": "PC-02",
            "check": "positive detector accepts linear current but not as kinetic",
            "computed": f"class={linear_current['classification']}; quad={len(linear_current['quad'])}",
            "expected": "linear_current with zero quadratic monomials",
            "status": "PASS" if linear_current["classification"] == "linear_current" and len(linear_current["quad"]) == 0 else "FAIL",
        },
        {
            "gate_id": "PC-03",
            "check": "native rows contain boundary-type linear currents",
            "computed": f"linear_rows={sum(1 for row in derivative_rows if row['classification'] == 'linear_current')}",
            "expected": "at least one linear current row; no spatial linear current",
            "status": "PASS" if any_linear and all_spatial_linear_zero else "FAIL",
        },
        {
            "gate_id": "NC-01",
            "check": "full 24-component torsion field theorem not constructed",
            "computed": nc_full_pde,
            "expected": "NOT_TRIGGERED",
            "status": nc_full_pde,
        },
    ]


def native_cell_ok(derivative_rows: list[dict[str, str]], topology: str, mode: str) -> tuple[bool, int]:
    rows = [row for row in derivative_rows if row["topology"] == topology and row["mode"] == mode]
    return bool(rows) and all(row["status"] == "PASS" for row in rows), len(rows)


def build_row(topology: str, mode: str, derivative_rows: list[dict[str, str]]) -> dict[str, str]:
    eta, v_var = sp.symbols("eta V", real=True)
    branch = derive_static_branch_function(topology, mode)
    max_degree = max_total_degree(branch, eta, v_var)
    shell = solve_auxiliary_shell(branch)
    solution = shell["solution"]
    aux_vars = shell["aux_vars"]
    shell_eta_zero = eta not in aux_vars or solution.get(eta, sp.S.Zero) == 0
    shell_v_zero = v_var not in aux_vars or solution.get(v_var, sp.S.Zero) == 0
    shell_ok = (
        shell["solution_status"] == "UNIQUE_NONDEGENERATE"
        and shell_eta_zero
        and shell_v_zero
    )
    native_ok, native_count = native_cell_ok(derivative_rows, topology, mode)
    status = "PASS" if max_degree <= 2 and shell_ok and native_ok else "FAIL"
    return {
        "topology": topology,
        "mode": mode,
        "max_degree_eta_V": str(max_degree),
        "shell_eta_zero": "n/a" if eta not in aux_vars else ("true" if shell_eta_zero else "false"),
        "shell_V_zero": "n/a" if v_var not in aux_vars else ("true" if shell_v_zero else "false"),
        "quadratic_derivative_absent": "true" if native_ok else "false",
        "native_derivative_rows": str(native_count),
        "native_gate_source": "in_process_bundled_native_ec_ny_engine",
        "status": status,
    }


def main() -> int:
    setup_log(__file__, log_dir=str(DATA_DIR))
    try:
        print("Paper07 torsion auxiliary field-promotion witness")
        print("=" * 65)
        print("full_24_component_torsion_claim=false")
        print("field_promotion_source=in-process bundled native EC+NY derivative computation")
        derivative_rows = native_derivative_rows()
        native_rows = detector_control_rows(derivative_rows)
        print_json_rows("native_derivative", derivative_rows, DERIVATIVE_FIELDS)
        print_json_rows("native_gate", native_rows, NATIVE_FIELDS)
        rows = []
        for topology in CANONICAL_TOPOLOGIES:
            for mode in MODES:
                if mode == "EH":
                    continue
                row = build_row(topology, mode, derivative_rows)
                rows.append(row)
                print(
                    f"  {topology}x{mode}: degree={row['max_degree_eta_V']} "
                    f"shell_eta={row['shell_eta_zero']} shell_V={row['shell_V_zero']} "
                    f"native_rows={row['native_derivative_rows']} "
                    f"status={row['status']}"
                )
        print_json_rows("torsion_cell", rows, ROW_FIELDS)
        native_gate_ok = all(row["status"] in {"PASS", "NOT_TRIGGERED"} for row in native_rows)
        print(f"native_gate_in_process_ok={native_gate_ok}")
        fail_count = sum(1 for row in rows if row["status"] != "PASS")
        fail_count += sum(1 for row in native_rows if row["status"] not in {"PASS", "NOT_TRIGGERED"})
        all_pass = fail_count == 0
        print()
        print(f"total_cells={len(rows)}")
        print(f"fail_count={fail_count}")
        if all_pass:
            print("overall_verdict=PASS_TORSION_AUXILIARY_FIELD_PROMOTION_WITNESS")
        else:
            print("overall_verdict=FAIL_TORSION_AUXILIARY_FIELD_PROMOTION_WITNESS")
        return 0 if all_pass else 1
    finally:
        teardown_log()


if __name__ == "__main__":
    raise SystemExit(main())

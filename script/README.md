# DPPUv7 paper07 Scripts Directory

-> [Japanese](README_ja.md)

**Paper**: "Dirac-before-Routh Constrained Spectra in Bianchi Class-A Einstein–Cartan–Nieh–Yan Cosmologies" (paper07)

Python packages and proof witness scripts for the paper07 Lorentzian EC+NY DBR constrained-spectrum analysis over the four Bianchi class-A backgrounds `T3`, `Nil3`, `Sol3`, and `S3`. The witnesses support rank-defect realization, scoped $k>0$ spectral completion, harmonic kernel correspondence, finite Palatini branch classification, torsion auxiliary/field-promotion checks, route comparison, and $S^3$ normalization invariance documented in the Japanese Markdown draft and Appendix B reproducibility notes.

Witness scripts use `script/dppu` .

---

## Directory Structure

```text
script/
|-- docs/                      # Technical documentation and conventions
|-- dppu/                      # Main Python package (DPPUv7)
|   |-- action/                # EC action, Lagrangian, reduced-sector helpers
|   |-- analysis/              # Dirac constraint analysis helpers
|   |-- curvature/             # Riemann, Ricci, spatial Lie-frame curvature
|   |-- dirac/                 # DBR/RBD constraint-closure engines
|   |-- engine/                # Metric, LC/EC connection, contortion, pipeline
|   |-- forms/                 # Exterior algebra and Nieh-Yan form
|   |-- kk/                    # Kaluza-Klein photon effective theory
|   |-- scanning/              # Parameter scan and diagnostics loader
|   |-- spectral/              # Harmonic bases, Hessians, Sol3 compact-edge spectrum
|   |-- topology/              # Class-A structure matrix and topology data
|   |-- torsion/               # Torsion modes, ansatz, Nieh-Yan density
|   `-- utils/                 # Common utilities and tee logging
|-- run_paper07.sh             # Batch runner for all paper07 witnesses
`-- scripts/                   # Execution scripts
    |-- paper07/               # paper07 proof witness scripts
    `-- visualize/             # Figure generation scripts
```

### `docs/` - Documentation

- [DPPUv7 Engine CONVENTIONS](docs/CONVENTIONS.md)
- [DPPUv7 SymPy guideline](docs/SymPy_guideline.md)
- [DPPUv7 LOGGING](docs/LOGGING.md)

---

## Package Overview (`dppu/`)

| Module | Role | Key Classes / Functions |
|--------|------|-------------------------|
| [`action`](dppu/action/README.md) | EC action, reduced sector, auxiliary shell | `reduced_sector`, `derive_static_branch_function`, `solve_auxiliary_shell` |
| [`analysis`](dppu/analysis/README.md) | Dirac–Bergmann and constraint helpers | `dirac` |
| [`curvature`](dppu/curvature/README.md) | Curvature tensors and spatial Lie-frame curvature | `RiemannTensor`, `spatial_lie` |
| [`dirac`](dppu/dirac/README.md) | DBR/RBD constraint closure and route comparison | `build_scalar_axis_projection`, `dirac_standard`, `routh_eliminate` |
| [`engine`](dppu/engine/README.md) | Metric, LC/EC connection, contortion, pipeline | `metric`, `levi_civita`, `ec_connection`, `contortion` |
| [`forms`](dppu/forms/README.md) | Exterior differential form algebra | `Form`, `wedge` |
| [`kk`](dppu/kk/README.md) | KK photon effective theory | `extractor` |
| [`scanning`](dppu/scanning/README.md) | Parameter scan and result loader | `parameter_scan`, `scan_results_loader` |
| [`spectral`](dppu/spectral/README.md) | Harmonic bases, reduced Hessians, Sol3 compact edge | `closure_matrix`, `scalar_axis_hessian_audits`, `compute_tables` |
| [`topology`](dppu/topology/README.md) | Class-A structure constants and invariants | `class_a`, `lz_invariants` |
| [`torsion`](dppu/torsion/README.md) | Torsion modes and Nieh-Yan variants | `Mode`, `construct_torsion_tensor` |
| [`utils`](dppu/utils/README.md) | Common utilities and tee logging | `setup_log`, `teardown_log`, `tee_logger` |

Shared witness helpers live in `scripts/paper07/paper07_common.py`.

---

## Script Overview (`scripts/`)

### `paper07/` - Proof Witness Scripts

| Script | Description |
|--------|-------------|
| `rank_defect_activation_witness.py` | Runs the `dppu.dirac` scalar-axis DBR closure engine and verifies $\mathrm{rank}(\mathcal A_{\rm DBR})=\mathrm{rank}(n)$, $d_{\rm zero}=3-\mathrm{rank}(n)$, and SO(3) frame rank invariance across four topologies and three scalar axes (Section 4, Appendix B). The formal $A=n$ identity is reported only as a notation check. |
| `spectral_completion_witness.py` | Checks H1–H3 positivity for $T^3$, $S^3$, $Nil^3$ (toroidal + Heisenberg towers), $Sol^3$ bulk AM–GM, and the `dppu.spectral.sol3_compact_edge` A-dual-orbit modified-Mathieu compact-edge $\lambda_1>0$ witness (Section 5, Appendix D). |
| `torsion_auxiliary_field_promotion_witness.py` | Verifies exact degree- $\le 2$ auxiliary shell and absence of quadratic derivative couplings under field promotion for AX/VT/MX across four topologies (Section 8, Appendix G). |
| `finite_fate_branch_witness.py` | Re-runs the `dppu.action.finite_reduction` Palatini symbolic series and checks exact-modulus / curvature-wall classification, tangent span, $T^3$ flat identity, and DeWitt trace-free norm (Section 7, Appendix F). |
| `harmonic_kernel_correspondence.py` | Uses `dppu.spectral` and `dppu.dirac` to verify closed 1-form space $=\ker n$ and the DBR-closure constraint-surface reduced Hessian kernel contact (Section 6, Appendix E). The old $H=n^Tn$ minimal model is not used. |
| `route_comparison_s3_branch_a_witness.py` | Recomputes DBR_KEEP_ALL and RBD Routh variant Nphys-by- $l$ / commutators from the same Branch-A local Lagrangian in `dppu.dirac` (Section 3, Appendix C). It also computes the $l=0$ RBD surviving scalar-pair $H_{PP}=-3/8$ directly from the Hamiltonian. Intermediate CSVs are not used as evidence. |
| `s3_normalization_factor_witness.py` | Under $C\to sC$, $\mathrm{rank}(n)$ and $\dim\ker n$ invariant; curvature and eigenvalue scale as $s^2$ (Section 2, Appendix A). |

### `dppu/` Engine Usage

The witnesses reconstruct the needed quantities locally from `script/dppu`: DBR/RBD route calculations from `dppu.dirac`, harmonic/reduced-Hessian/Sol3 compact-edge calculations from `dppu.spectral`, and finite Palatini branch plus torsion auxiliary checks from `dppu.action`.

#### `paper07/` Execution Matrix

Each witness writes to `../data/<script_name>_YYYYMMDD_HHMMSS.log` when run from `script/`. Set `DPPU_LOG_DIR` to override the log directory (`run_paper07.sh` sets this automatically).

| Script | Command from `script/` | Expected verdict |
|--------|------------------------|------------------|
| `rank_defect_activation_witness.py` | `python scripts/paper07/rank_defect_activation_witness.py` | `overall_verdict=PASS_RANK_DEFECT_ACTIVATION_WITNESS` |
| `spectral_completion_witness.py` | `python scripts/paper07/spectral_completion_witness.py` | `overall_verdict=PASS_SPECTRAL_COMPLETION_WITNESS` |
| `torsion_auxiliary_field_promotion_witness.py` | `python scripts/paper07/torsion_auxiliary_field_promotion_witness.py` | `overall_verdict=PASS_TORSION_AUXILIARY_FIELD_PROMOTION_WITNESS` |
| `finite_fate_branch_witness.py` | `python scripts/paper07/finite_fate_branch_witness.py` | `overall_verdict=PASS_FINITE_FATE_BRANCH_WITNESS` |
| `harmonic_kernel_correspondence.py` | `python scripts/paper07/harmonic_kernel_correspondence.py` | `overall_verdict=PASS_HARMONIC_KERNEL_CORRESPONDENCE` |
| `route_comparison_s3_branch_a_witness.py` | `python scripts/paper07/route_comparison_s3_branch_a_witness.py` | `overall_verdict=PASS_ROUTE_COMPARISON_S3_BRANCH_A_WITNESS` |
| `s3_normalization_factor_witness.py` | `python scripts/paper07/s3_normalization_factor_witness.py` | `overall_verdict=PASS_S3_NORMALIZATION_FACTOR_WITNESS` |

Batch run (proposal order):

```bash
bash run_paper07.sh
# optional: bash run_paper07.sh --output-dir /path/to/data
```

---

### `visualize/` - Figures

Figure scripts do not create execution logs. Each script writes to a fixed path under `../LaTeX/figures/`.

| File | Description | Output |
|------|-------------|--------|
| `fig04_spectral_completion_map.py` | Four-panel scoped $k>0$ spectral completion map: $T^3$ Fourier, $S^3$ TT harmonics, $Nil^3$ toroidal + Heisenberg tower, $Sol^3$ bulk/edge split (Figure 4, Section 5) | `../LaTeX/figures/fig04_spectral_completion_map.svg` |

---

## Quick Start

```bash
# From the script directory
cd script

pip install -r requirements.txt

# Run all witnesses
bash run_paper07.sh

# Or run one witness
python scripts/paper07/rank_defect_activation_witness.py

# Regenerate Figure 4
python scripts/visualize/fig04_spectral_completion_map.py
```

PowerShell:

```powershell
Set-Location script
pip install -r requirements.txt
bash run_paper07.sh
python scripts/visualize/fig04_spectral_completion_map.py
```

---

## Paper07 Reproducibility Map

| Paper location | Supporting script |
|----------------|-------------------|
| Section 2, Appendix A: $S^3$ normalization and scale invariance | `scripts/paper07/s3_normalization_factor_witness.py` |
| Section 3, Appendix C: DBR vs RBD negative control | `scripts/paper07/route_comparison_s3_branch_a_witness.py` |
| Section 4, Appendix B: rank-defect activation and frame invariance | `scripts/paper07/rank_defect_activation_witness.py` |
| Section 5, Appendix D: H1/H2/H3 spectral completion | `scripts/paper07/spectral_completion_witness.py` |
| Section 6, Appendix E: harmonic 1-form / Hessian kernel | `scripts/paper07/harmonic_kernel_correspondence.py` |
| Section 7, Appendix F: finite modulus/wall split | `scripts/paper07/finite_fate_branch_witness.py` |
| Section 8, Appendix G: torsion auxiliary and field promotion | `scripts/paper07/torsion_auxiliary_field_promotion_witness.py` |
| Appendix B §B.6: batch reproducibility | `run_paper07.sh` + all `scripts/paper07/*_witness.py` logs |

## Scope Notes

- All physical claims in the draft are under the **DBR (Dirac-before-Routh)** protocol unless explicitly marked as route comparison.
- H1–H3 spectral results are scoped to **analyzed harmonic sectors**; they are not full field-theory PDE spectral theorems.
- $Sol^3$ compact-edge positivity is a **representative numeric** witness, not an all-Anosov closed-form theorem.
- Finite modulus/wall classification is scoped to **finite Palatini branch / homogeneous minisuperspace**; finite moduli count is generally not equal to $d_{\rm zero}$.
- Torsion robustness applies to **NY-MAIN route** and **AX/VT/MX ansatz scalars** only; it is not a full 24-component torsion field-theory theorem.
- RBD route discrepancy is a **route-level artifact**, not a claim that RBD always fails or that an extra scalar is a physical ghost.
- Witness scripts reconstruct definitions from the paper and `dppu` locally.

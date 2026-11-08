# Paper07 scientific supplement — v2

The package derives the specified formulas, coefficients and finite-time examples. The paper gives the hypotheses, constructions and proofs; passing finite checks alone is not a theorem of nonlinear existence or general stability.

## Commands

From the repository root, Python 3.11+:

```sh
python -m pip install -r script/requirements.txt
python -B script/paper07.py recompute --output results_base_new
python -B script/perturbations.py recompute --model all --output results_models_new
python -B script/perturbations.py export --input results_models_new --output selected_data_new
python -B script/paper07.py render --data data/revised --output figures_base_new
python -B script/perturbations.py render --data data/perturbations --output figures_models_new
```

Every output directory must be new. `recompute` generates scientific results, `export` selects publication data, and `render` reads saved numeric arrays without importing symbolic model producers. Both renderers verify data hashes; together they produce the five PNG/PDF figure pairs. Replacing `data/perturbations` with `selected_data_new` renders a fresh computation.

Models can be run separately:

```sh
python -B script/perturbations.py recompute --model sol --output sol_new
python -B script/perturbations.py recompute --model lorenz --canonical-input sol_new/reduction_core.json --output lorenz_new
```

With `--model all`, the Lorenz calculation uses the freshly computed Sol canonical map. A standalone Lorenz run without `--canonical-input` reads the shipped map. The fixed `data/initial/sol_profile.json` is checked against its manifest and used unchanged. For standalone model outputs, collect the six named subdirectories before using `export`.

The model driver sets the SymPy Python ground backend, one BLAS/OpenMP thread and a fixed Python hash seed for its subprocesses. The original formula command also supports setting `SYMPY_GROUND_TYPES=python` explicitly. No other repository, work folder, network access, credentials or external audit tools are needed.

The complete model run generates roughly 0.6 GB, primarily four Sol trajectories and their quadrature fields. These can require substantially more RAM than their file size while JSON is parsed. Runtime depends on symbolic backend and hardware; progress is printed per producer. The distributed selected data are much smaller. Intermediate Nil oscillator buffers and full Sol trajectories are generated on demand, not included as reader-facing saved artifacts. Do not interrupt a run and reuse its output path: keep the partial directory and choose a new one.

## Models and authority

| Model | Definitions and outputs | Manuscript |
| --- | --- | --- |
| `paper07.py` | Full contortion form, homogeneous regular constraints, scalar controls, torus quadratic action, geometric comparisons | Sections 2–5; Appendices A–E and I |
| `flat` | Four KIDs and AMM cone normalization; explicit single/two-wave corrections; 40/34/10 configuration connection comparison | Sections 4, 6–7; Appendices B–D |
| `backgrounds` | Exact Nil/Sol/Taub vacua and quotients; Nil base and Sol Gowdy reductions, endpoint controls, examples | Appendices F–H |
| `nil` | Central-nonzero tensor/oscillator closure; five explicit grades, L=2 Einstein reconstruction and example | Appendices F–H |
| `sol` | Infinite-orbit differential operators, right inverse, canonical embedding, independent Einstein and Cartan checks | Appendices F–G |
| `taub` | Four low Wigner blocks; j=3/2 reconstruction, Cartan/endpoints, KIDs and example | Appendices F–H |
| `lorenz` | Expanded reduced fourth-order operator (SG-05), independent Lorenz system, common initial maps, four trajectories, comparisons and actual dense-polynomial residuals | Appendix H |

The underlying `paper07_results` modules use signature (−,+,+,+), epsilon_0123=+1 and geometric torsion T=de+omega wedge e. The structure equation, fixed quotient, physical spherical radius and doubled off-diagonal canonical momenta are in Appendix A. Auxiliary connection comparisons use the same complete NY current and corresponding endpoint representatives.

`quadratic_dirac` accepts constant-coefficient quadratic finite-dimensional actions; it is not a nonlinear Dirac solver. The curved producers build the full linear ADM constraints from the displayed background and tensor basis. Formal differential adjoints in the Sol reduction use the common smooth core in Appendix F, not a claimed closed self-adjoint extension.

## Reading the data

Symbolic expressions use SymPy strings; matrices are row-major nested lists. Sol differential operators are lists of coefficient matrices multiplying successive derivatives on the **right**: sum A_j(t,xi) d_xi^j. Composition and adjoints differentiate coefficients. The coefficient symbols F, K, M represent exp(t²)/sqrt(t), k0 exp(−xi), l0 exp(xi). Saved names keep the code symbols: `exp_2P0` in `backgrounds/nil_base.json` is ζ = exp(2P̄) of (NG-02), and `L` in `backgrounds/sol_gowdy_example.json` is the Sol period L_Sol. The real/complex phases and component ordering are part of each result and Appendices A, G–H.

Model JSON files retain exact algebraic checks, coefficients and scope. The finite checks support the paper's derivations; the flat-torus AMM sufficiency claim and scoped Sol no-KID sufficiency claim come from the cited theorems and their function-space arguments. The Sol result includes admissible real reconstructed initial data, not arbitrary reduced L2 inputs or a well-posedness theorem for SG-05. Read [the scope table](../data/perturbations/SCOPE.md) with the outputs. In particular, a small four-dimensional Einstein component residual does not certify convergence of the reduced fourth-order system or a whole-line solution.

`dense_residual.py` differentiates SciPy's actual seventh-degree DOP853 interpolant and stores the implementation source and SciPy version used. It does not replace time derivatives by the ODE right-hand side. This use of SciPy's private dense-output class was tested with SciPy 1.16 (`requirements.txt` bounds it to 1.16.x); review compatibility before using a later SciPy.

The v1 `dppu` API and output schemas are superseded. Formula data are in `data/revised`; the perturbation models have a separate manifest in `data/perturbations`. See [the comparison guide](legacy_comparison.md).

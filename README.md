# Lorentzian Einstein–Cartan Reduction on Compact Homogeneous Geometries: Constraints and Inhomogeneous Perturbations

Muacca · **Version 2** · Revised 2026-11-08 · First release 2026-06-14

We study inhomogeneous constraints in source-free Lorentzian Einstein–Cartan theory with constant complete Nieh–Yan coupling on fixed compact geometries. The regular homogeneous metric counts remain 5, 3, 2 and 2 pairs on the torus, Nil, Sol and the three-sphere, with symmetric-background exceptions treated separately. For nonzero flat-torus modes, connection-retaining and metric reductions agree in the action and reduced two-form. Under the stated compact CMC hypotheses, the Arms–Marsden–Moncrief theorem makes the four Taub conditions necessary and sufficient for local integrability; explicit standing-wave and nonparallel two-wave corrections satisfy the local vacuum equations through second order. On selected time-dependent curved vacuum backgrounds, we construct closed tensor blocks on LRS Nil, infinite dual-orbit systems on Sol, and low Wigner blocks on Taub, with gauges, canonical reductions and reconstruction. Background-dependent quadratic KID conditions are distinguished from nonlinear lifting. Finite-time examples include an independent Lorenz evolution on one Sol orbit, with interior spatial, window and time convergence comparisons. Convergence in the reduced fourth-order form, whole-line error bounds, general stability and curved-background lifting sufficiency are not established.

## Read

- [English PDF](DPPUv7-paper07.pdf) and [LaTeX source](LaTeX/main.tex)
- [日本語本文・目次](README_ja.md)
- [Revision history](CHANGELOG.md) and [references](REFERENCES.md)
- [Code and commands](script/README.md), [data guide](data/README.md), [adopted and unresolved scope](data/perturbations/SCOPE.md)

The regular homogeneous counts, field species and complete real-block totals use different units. Nil's lower-grade missing components differ from Taub's KID rank exceptions. The flat-torus integrability theorem, explicit second-order constructions, curved linear reductions and numerical examples are separate levels of evidence. Geometric Laplacians in Appendix E are not a common physical spectrum.

## Reproduce

From the repository root, with Python 3.11 or later:

```sh
python -m pip install -r script/requirements.txt
python -B script/paper07.py recompute --output results_base_new
python -B script/perturbations.py recompute --model all --output results_models_new
python -B script/paper07.py render --data data/revised --output figures_base_new
python -B script/perturbations.py render --data data/perturbations --output figures_models_new
```

All output directories must be new. `all` includes the four large Sol trajectories and their dense-time residual check; models can also run separately. Rendering reads saved data and does not rerun derivations. See the script guide for dependencies, output sizes and validation scope. Build the PDF by running `pdflatex -interaction=nonstopmode -halt-on-error main.tex` twice inside `LaTeX`; the supplied figures already match the saved data.

Use the matching v2 manuscript, code and data. This supplement supersedes the v1 local engine and its saved outputs; the [comparison guide](script/legacy_comparison.md) relates the v1 conventions to the present ones. Correction and withdrawal notes are in Appendix I and the changelog.

## Availability and acknowledgment

[Repository](https://github.com/Muacca/DPPUv7-paper07) · [Existing Zenodo record](https://doi.org/10.5281/zenodo.20678646). Cite the edition and date of the artifact actually obtained. [MIT License](LICENSE.txt). Contact: [muacca@dmwp.jp](mailto:muacca@dmwp.jp).

The original manuscript acknowledged Claude Opus 4.8, Gemini 3.1 Pro, Grok 4.3 and GPT-5.5 Thinking (June 2026). OpenAI Codex assisted this revision with symbolic derivations, code, Japanese–English exposition and deterministic plotting. The author is responsible for the content and remaining errors.

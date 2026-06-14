# Dirac-before-Routh Constrained Spectra in Bianchi Class-A Einstein–Cartan–Nieh–Yan Cosmologies: Rank-Defect Kernels, Spectral Robustness, and Torsion-Sector Closure

[日本語 README](README_ja.md)

## Paper Contents

- [Official PDF (English, built from LaTeX)](DPPUv7-paper07.pdf)
- [Draft in Japanese Markdown](DPPUv7-paper07_sec01.md)

## Directory Structure

- `DPPUv7-paper07_sec01.md` … `DPPUv7-paper07_sec09.md` — Main sections (Japanese Markdown draft)
- `DPPUv7-paper07_appA.md` … `DPPUv7-paper07_appG.md` — Appendices
- `LaTeX/` — LaTeX assets
  - `figures/` — Figures (e.g. Figure 4 spectral completion map)
- `script/` — Proof witness scripts and the `dppu` Python package
  - For script details, see [script/README.md](script/README.md)
- `data/` — Timestamped witness logs (`<script_name>_YYYYMMDD_HHMMSS.log`)

## Building the LaTeX PDF

When `LaTeX/main.tex` is added, build from the LaTeX directory with **at least two** `pdflatex` passes so cross-references resolve correctly.

```bash
cd LaTeX
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

### Troubleshooting

Check the log on error:

```bash
grep -i "error\|warning" main.log | tail -20
```

Clean rebuild:

```bash
cd LaTeX
rm -f *.aux *.log *.out *.synctex.gz *.fdb_latexmk *.fls
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

## Reproducibility (Witness Scripts)

From `script/`:

```bash
pip install -r requirements.txt
bash run_paper07.sh
```

Logs are written to `../data/`. Witness scripts do not read external analysis-phase CSV/JSON; they reconstruct paper definitions locally. See [script/README.md](script/README.md) for the paper-to-script map.

## Author

Muacca

## License

See the LICENSE file in the repository root.
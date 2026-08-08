# Dirac-before-Routh Constrained Spectra in Bianchi Class-A Einstein–Cartan–Nieh–Yan Cosmologies: Rank-Defect Kernels, Spectral Robustness, and Torsion-Sector Closure

[10.5281/zenodo.20678646](https://doi.org/10.5281/zenodo.20678646)

[English README](README.md)

## 論文内容

- [正式な PDF（英語。LaTeX からビルドしたもの）](DPPUv7-paper07.pdf)
- [ドラフトの日本語版 Markdown](DPPUv7-paper07_sec01.md)


## ディレクトリ構造

- `DPPUv7-paper07_sec01.md` … `DPPUv7-paper07_sec09.md` — 本文（日本語 Markdown ドラフト）
- `DPPUv7-paper07_appA.md` … `DPPUv7-paper07_appG.md` — 付録
- `LaTeX/` — LaTeX 関連ファイル
  - `figures/` — 図（例: Figure 4 spectral completion map）
- `script/` — proof witness スクリプトと `dppu` Python パッケージ
  - スクリプトの詳細は [script/README_ja.md](script/README_ja.md) を参照ください。
- `data/` — witness 実行ログ（`<script_name>_YYYYMMDD_HHMMSS.log`）

## LaTeX のビルド方法

`LaTeX/main.tex` が追加されたら、LaTeX ディレクトリで pdflatex を実行します。相互参照（`\ref`、`\label` など）を正しく解決するため、**最低 2 回のコンパイルが必要**です。

```powershell
cd LaTeX
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

### トラブルシューティング

エラー時はログを確認：

```powershell
cat main.log | Select-String -Pattern "Error|Warning" | Select-Object -Last 20
```

クリーンビルド：

```powershell
cd LaTeX
Remove-Item *.aux, *.log, *.out, *.synctex.gz, *.fdb_latexmk, *.fls -ErrorAction SilentlyContinue
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

## 再現性（witness スクリプト）

`script/` から：

```powershell
pip install -r requirements.txt
bash run_paper07.sh
```

ログは `../data/` に出力されます。witness スクリプトは外部の解析フェーズ CSV/JSON を読まず、論文内の定義をローカルに再構成します。論文との対応表は [script/README_ja.md](script/README_ja.md) を参照してください。

## 著者

Muacca

## ライセンス

リポジトリ直下の LICENSE ファイルを参照してください。
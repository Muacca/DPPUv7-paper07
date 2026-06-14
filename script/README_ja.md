# DPPUv7 paper07 スクリプトディレクトリ

-> [English](README.md)

**論文**: "Dirac-before-Routh Constrained Spectra in Bianchi Class-A Einstein–Cartan–Nieh–Yan Cosmologies"（paper07）

Bianchi class-A 背景 `T3`, `Nil3`, `Sol3`, `S3` 上の Lorentzian EC+NY DBR 制約スペクトル解析のための Python パッケージ群と proof witness スクリプトです。rank-defect realization、scoped $k>0$ spectral completion、harmonic kernel 対応、finite Palatini branch 分類、torsion auxiliary/field-promotion 検査、route 比較、 $S^3$ 正規化不変性を、日本語 Markdown ドラフトおよび Appendix B の再現性メモと対応づけます。

witness スクリプトは `script/dppu` を利用します。

---

## ディレクトリ構成

```text
script/
|-- docs/                      # 技術ドキュメントと規約
|-- dppu/                      # メイン Python パッケージ（DPPUv7）
|   |-- action/                # EC action, Lagrangian, reduced-sector helper
|   |-- analysis/              # Dirac 制約解析 helper
|   |-- curvature/             # Riemann, Ricci, spatial Lie-frame curvature
|   |-- dirac/                 # DBR/RBD constraint-closure engine
|   |-- engine/                # metric, LC/EC connection, contortion, pipeline
|   |-- forms/                 # exterior algebra と Nieh-Yan form
|   |-- kk/                    # Kaluza-Klein photon effective theory
|   |-- scanning/              # parameter scan と diagnostics loader
|   |-- spectral/              # harmonic basis, Hessian, Sol3 compact-edge spectrum
|   |-- topology/              # class-A 構造行列と topology data
|   |-- torsion/               # torsion mode, ansatz, Nieh-Yan density
|   `-- utils/                 # 共通 utility と tee logging
|-- run_paper07.sh             # paper07 witness 一括実行
`-- scripts/                   # 実行スクリプト
    |-- paper07/               # paper07 proof witness スクリプト
    `-- visualize/             # 図表生成スクリプト
```

### `docs/` - ドキュメント

- [DPPUv7 Engine CONVENTIONS](docs/CONVENTIONS_ja.md)
- [DPPUv7 SymPy guideline](docs/SymPy_guideline_ja.md)
- [DPPUv7 LOGGING](docs/LOGGING_ja.md)

---

## パッケージ概要（`dppu/`）

| モジュール | 役割 | 主なクラス / 関数 |
|------------|------|-------------------|
| [`action`](dppu/action/README_ja.md) | EC action, reduced sector, auxiliary shell | `reduced_sector`, `derive_static_branch_function`, `solve_auxiliary_shell` |
| [`analysis`](dppu/analysis/README_ja.md) | Dirac–Bergmann と制約 helper | `dirac` |
| [`curvature`](dppu/curvature/README_ja.md) | curvature tensor と spatial Lie-frame curvature | `RiemannTensor`, `spatial_lie` |
| [`dirac`](dppu/dirac/README_ja.md) | DBR/RBD constraint closure と route comparison | `build_scalar_axis_projection`, `dirac_standard`, `routh_eliminate` |
| [`engine`](dppu/engine/README_ja.md) | metric, LC/EC connection, contortion, pipeline | `metric`, `levi_civita`, `ec_connection`, `contortion` |
| [`forms`](dppu/forms/README_ja.md) | exterior differential form algebra | `Form`, `wedge` |
| [`kk`](dppu/kk/README_ja.md) | KK photon effective theory | `extractor` |
| [`scanning`](dppu/scanning/README_ja.md) | parameter scan と result loader | `parameter_scan`, `scan_results_loader` |
| [`spectral`](dppu/spectral/README_ja.md) | harmonic basis, reduced Hessian, Sol3 compact edge | `closure_matrix`, `scalar_axis_hessian_audits`, `compute_tables` |
| [`topology`](dppu/topology/README_ja.md) | class-A 構造定数と不変量 | `class_a`, `lz_invariants` |
| [`torsion`](dppu/torsion/README_ja.md) | torsion mode と Nieh-Yan variant | `Mode`, `construct_torsion_tensor` |
| [`utils`](dppu/utils/README_ja.md) | 共通 utility と tee logging | `setup_log`, `teardown_log`, `tee_logger` |

witness 共通 helper は `scripts/paper07/paper07_common.py` にあります。

---

## 実行スクリプト概要（`scripts/`）

### `paper07/` - proof witness スクリプト

| スクリプト | 説明 |
|------------|------|
| `rank_defect_activation_witness.py` | `dppu.dirac` の scalar-axis DBR closure engine を実行し、四 topology × 三 scalar axis で $\mathrm{rank}(\mathcal A_{\rm DBR})=\mathrm{rank}(n)$、 $d_{\rm zero}=3-\mathrm{rank}(n)$、SO(3) frame における rank 不変を確認（Section 4, Appendix B）。形式的な $A=n$ 恒等式は記法チェックとしてのみ併記する。 |
| `spectral_completion_witness.py` | $T^3$/ $S^3$/ $Nil^3$（toroidal + Heisenberg）/ $Sol^3$ bulk AM–GM と、`dppu.spectral.sol3_compact_edge` の A-dual-orbit modified-Mathieu compact-edge $\lambda_1>0$ witness を確認（Section 5, Appendix D）。 |
| `torsion_auxiliary_field_promotion_witness.py` | AX/VT/MX × 四 topology で exact degree $\le 2$ auxiliary shell と二次微分結合不在を確認（Section 8, Appendix G）。 |
| `finite_fate_branch_witness.py` | `dppu.action.finite_reduction` の Palatini symbolic series を再実行し、finite Palatini branch 上の exact modulus / curvature wall 分類、tangent span、 $T^3$ flat 恒等式、DeWitt trace-free norm を確認（Section 7, Appendix F）。 |
| `harmonic_kernel_correspondence.py` | `dppu.spectral` と `dppu.dirac` を使い、閉 1-form 空間 $=\ker n$ と DBR closure 後の constraint-surface reduced Hessian kernel 接触を確認（Section 6, Appendix E）。旧 $H=n^Tn$ minimal model は使わない。 |
| `route_comparison_s3_branch_a_witness.py` | `dppu.dirac` の同一 Branch-A local Lagrangian から DBR_KEEP_ALL と RBD Routh variants の Nphys-by- $l$ と commutator を再計算する（Section 3, Appendix C）。 $l=0$ RBD surviving scalar pair の $H_{PP}=-3/8$ も Hamiltonian から直接計算する。中間 CSV は証拠として使わない。 |
| `s3_normalization_factor_witness.py` | $C\to sC$ で $\mathrm{rank}(n)$/ $\dim\ker n$ 不変、曲率・固有値は $s^2$ スケールを確認（Section 2, Appendix A）。 |

### `dppu/` engine 利用方針

witness は `script/dppu` の各 engine からローカルに再構成します。DBR/RBD route 計算は `dppu.dirac`、harmonic/reduced-Hessian/Sol3 compact-edge 計算は `dppu.spectral`、finite Palatini branch と torsion auxiliary checks は `dppu.action` を入口にします。

#### `paper07/` 実行仕様

各 witness は `script/` から実行した場合、`../data/<script_name>_YYYYMMDD_HHMMSS.log` に log を出力します。`run_paper07.sh` が `DPPU_LOG_DIR` を自動設定します。

| スクリプト | `script/` からの実行コマンド | 期待される verdict |
|------------|------------------------------|--------------------|
| `rank_defect_activation_witness.py` | `python scripts/paper07/rank_defect_activation_witness.py` | `overall_verdict=PASS_RANK_DEFECT_ACTIVATION_WITNESS` |
| `spectral_completion_witness.py` | `python scripts/paper07/spectral_completion_witness.py` | `overall_verdict=PASS_SPECTRAL_COMPLETION_WITNESS` |
| `torsion_auxiliary_field_promotion_witness.py` | `python scripts/paper07/torsion_auxiliary_field_promotion_witness.py` | `overall_verdict=PASS_TORSION_AUXILIARY_FIELD_PROMOTION_WITNESS` |
| `finite_fate_branch_witness.py` | `python scripts/paper07/finite_fate_branch_witness.py` | `overall_verdict=PASS_FINITE_FATE_BRANCH_WITNESS` |
| `harmonic_kernel_correspondence.py` | `python scripts/paper07/harmonic_kernel_correspondence.py` | `overall_verdict=PASS_HARMONIC_KERNEL_CORRESPONDENCE` |
| `route_comparison_s3_branch_a_witness.py` | `python scripts/paper07/route_comparison_s3_branch_a_witness.py` | `overall_verdict=PASS_ROUTE_COMPARISON_S3_BRANCH_A_WITNESS` |
| `s3_normalization_factor_witness.py` | `python scripts/paper07/s3_normalization_factor_witness.py` | `overall_verdict=PASS_S3_NORMALIZATION_FACTOR_WITNESS` |

一括実行（提案順）:

```bash
bash run_paper07.sh
# 任意: bash run_paper07.sh --output-dir /path/to/data
```

---

### `visualize/` - 図表

図生成スクリプトは実行ログを作成しません。各スクリプトは `../LaTeX/figures/` 配下の固定パスに出力します。

| ファイル | 説明 | 出力先 |
|----------|------|--------|
| `fig04_spectral_completion_map.py` | scoped $k>0$ spectral completion の 4 パネル模式図： $T^3$ Fourier、 $S^3$ TT harmonics、 $Nil^3$ toroidal + Heisenberg tower、 $Sol^3$ bulk/edge 分離（Figure 4, Section 5） | `../LaTeX/figures/fig04_spectral_completion_map.svg` |

---

## クイックスタート

```bash
# script ディレクトリへ移動
cd script

pip install -r requirements.txt

# 全 witness を実行
bash run_paper07.sh

# 単体実行
python scripts/paper07/rank_defect_activation_witness.py

# Figure 4 を再生成
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

## paper07 再現性マップ

| 論文内の位置 | 対応スクリプト |
|--------------|----------------|
| Section 2, Appendix A: $S^3$ 正規化とスケール不変性 | `scripts/paper07/s3_normalization_factor_witness.py` |
| Section 3, Appendix C: DBR vs RBD negative control | `scripts/paper07/route_comparison_s3_branch_a_witness.py` |
| Section 4, Appendix B: rank-defect activation と frame 不変性 | `scripts/paper07/rank_defect_activation_witness.py` |
| Section 5, Appendix D: H1/H2/H3 spectral completion | `scripts/paper07/spectral_completion_witness.py` |
| Section 6, Appendix E: harmonic 1-form / Hessian kernel | `scripts/paper07/harmonic_kernel_correspondence.py` |
| Section 7, Appendix F: finite modulus/wall split | `scripts/paper07/finite_fate_branch_witness.py` |
| Section 8, Appendix G: torsion auxiliary と field promotion | `scripts/paper07/torsion_auxiliary_field_promotion_witness.py` |
| Appendix B §B.6: 一括再現性 | `run_paper07.sh` + 全 `scripts/paper07/*_witness.py` log |

## スコープ注記

- ドラフトの物理主張は、route 比較を除き **DBR（Dirac-before-Routh）プロトコル**下の結果です。
- H1–H3 は **解析済み harmonic sector** に限定され、full field-theory PDE spectral theorem ではありません。
- $Sol^3$ compact edge の正性は **代表的数値** witness であり、all-Anosov 閉形式定理ではありません。
- finite modulus/wall 分類は **finite Palatini branch / homogeneous minisuperspace** 内に限定され、finite moduli 数は一般に $d_{\rm zero}$ と一致しません。
- torsion 頑健性は **NY-MAIN route** と **AX/VT/MX ansatz scalars** に限定され、full 24 成分 torsion 場理論の定理ではありません。
- RBD の余分 scalar は **route-level artifact** であり、「RBD は常に失敗する」や physical ghost の主張ではありません。
- witness は、論文定義と `dppu` からローカルに再構成します。

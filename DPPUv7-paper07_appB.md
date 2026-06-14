# DPPUv7-paper07 — Appendix B: Native Local Action and DBR Computational Checks

[← Appendix A](DPPUv7-paper07_appA.md) | [Appendix C →](DPPUv7-paper07_appC.md)

---

## Appendix B. Native local action and DBR computational checks

### B.1 Reduced action ingredients

homogeneous 背景 ansatz（[Section 2.3](DPPUv7-paper07_sec02.md)）の下で、first-order EC+NY 作用は有限自由度の reduced Lagrangian に落ちる。独立変数は scale factors $q_a(t)$、lapse $N(t)$、homogeneous 接続成分、torsion ansatz scalars $(\eta, V)$ である。摂動解析では、この reduced 系を scalar / spin-1 / spin-2 sector にブロック分解した二次作用を用いる。

### B.2 Symbolic Dirac–Bergmann engine

制約解析は symbolic 計算エンジンで実行した。手続きは標準の Dirac–Bergmann アルゴリズム [10–12] に忠実である:

1. Legendre 変換と primary constraints の抽出。
2. consistency 条件による secondary constraints の生成（閉包まで反復）。
3. first/second class の分類と、second class 対への Dirac bracket 構成。
4. 閉包後の Routh 型消去による reduced スペクトルの提示（DBR の順序）。

エンジンの判定は特定の topology の期待値を事前に固定しない。実装上は各 topology の local builder から scalar-axis DBR 閉包を走らせ、primary preservation、secondary generation、constraint Poisson matrix rank、first/second class split から activation entry を読む。構造定数と構造行列の不定記号レベルの関係は formal identity として別に確認し、DBR activation としての物理的意味は axis-by-axis closure audit で確認する。

**§4.3 proof skeleton との対応.** scalar-axis 摂動 $\boldsymbol\beta\in\mathbb R^3$ に対する secondary activation 写像を $\mathcal A_{\rm DBR}$ と書く。class-A Milnor frame では、構造写像の formal identity として

$$
\mathcal A_{\rm form} = n
$$

が成り立つ。この恒等式は記法・構造写像の確認であり、DBR closure を代替しない。`rank_defect_activation_witness.py` は、同梱の DBR 制約閉包エンジン（`dppu.dirac`）を呼び、四 topology × 三 scalar axis について $\mathcal A_{\rm DBR}$ の entry を DBR closure から計算する。その結果、engine-derived activation matrix は $\mathrm{rank}(\mathcal A_{\rm DBR})=\mathrm{rank}(n)$ を満たし、 $\dim\ker(\mathcal A_{\rm DBR})=d_{\rm zero}$ を与える。

### B.3 Rank verification conventions

制約行列・kinetic Hessian の rank 判定には以下の冗長化を課す:

- **二重経路**: symbolic rank と数値 rank の一致を要求する。
- **二点 generic 検査**: 数値 rank は 2 つの generic な背景点で評価し、一致を要求する（単一点の偶然退化を排除）。
- **mod- $p$ 検査**: 大型 symbolic 行列には素数 modulus での rank 計算を併用し、桁落ちによる偽 rank を排除する。

### B.4 Frame-invariance check protocol

[Section 4.4](DPPUv7-paper07_sec04.md) の frame 不変性は、まず $n\mapsto RnR^{\mathsf T}$ の合同変換で $\mathrm{rank}(n)$ と $\dim\ker n$ が保存される解析的議論に基づく。数値検証は次の二重経路で補強した: (a) ランダム SO(3) 回転（20 回転すべてで $\dim\ker n$ 保存）、(b) 回転後の構造定数（テンソル変換）からの $n$ の再構成と kernel 次元の再計算。`rank_defect_activation_witness.py` でも有理直交行列 $R$ による $\mathrm{rank}(RnR^{\mathsf T})=\mathrm{rank}(n)$ を確認する。

### B.5 同一プロトコル下の比較手続き

四 topology の比較はすべて同一のエンジン・同一の sector 分解・同一の rank 検証規約の下で行った。 $S^3$ は手法の cross-check を兼ね、scalar-axis 射影の null 検査（等方 $S^3$ で全軸 0）を通過している。

### B.6 Reproducibility notes

論文内の主要式・判定は `script/scripts/paper07/` 配下の独立 witness スクリプトで再現可能である。

| witness | 対応節 |
|---|---|
| `rank_defect_activation_witness.py` | §4.3–4.4, B.2–B.4（formal identity は notation check、主証拠は `dppu.dirac` の DBR closure audit） |
| `spectral_completion_witness.py` | §5, D |
| `torsion_auxiliary_field_promotion_witness.py` | §8, G |
| `finite_fate_branch_witness.py` | §7, F |
| `harmonic_kernel_correspondence.py` | §6, E |
| `route_comparison_s3_branch_a_witness.py` | §3.3, C |
| `s3_normalization_factor_witness.py` | §2.2, A.3 |

---

[← Appendix A](DPPUv7-paper07_appA.md) | [Appendix C →](DPPUv7-paper07_appC.md)

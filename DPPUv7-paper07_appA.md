# DPPUv7-paper07 — Appendix A: Conventions

[← Section 9](DPPUv7-paper07_sec09.md) | [Appendix B →](DPPUv7-paper07_appB.md)

---

## Appendix A. Conventions

### A.1 Signature, Hodge, and $\epsilon$ conventions

- metric signature: $(-,+,+,+)$（Lorentzian）。
- internal indices $I,J,\dots\in\{0,1,2,3\}$、空間 frame indices $a,b,\dots\in\{1,2,3\}$。
- Levi-Civita 記号: $\epsilon_{0123}=+1$。空間 3 次元では $\epsilon_{123}=+1$。
- Hodge dual は Lorentzian 規約で定義し、4-form の向きは $e^0\wedge e^1\wedge e^2\wedge e^3$ を正とする。
- curvature/torsion: $R^{IJ}=d\omega^{IJ}+\omega^I{}_K\wedge\omega^{KJ}$, $T^I=de^I+\omega^I{}_J\wedge e^J$。

### A.2 Structure matrix sign convention

構造方程式 $d\sigma^i=-\tfrac12 C^i{}_{jk}\sigma^j\wedge\sigma^k$ に対し、

$$
n^{ij}=\tfrac12\,\epsilon^{ikl}\,C^{j}{}_{kl}
$$

と定義する（class A で対称）。Milnor frame での対角値の符号は Table 1 の規約（ $Nil^3$: $\mathrm{diag}(1,0,0)$、 $Sol^3$: $\mathrm{diag}(1,-1,0)$、 $S^3$: $\mathrm{diag}(1,1,1)$ ）に固定する。 $\mathrm{rank}(n)$ と $\dim\ker n$ は対角値の正規化（スケール）に依存しない。

### A.3 $S^3$ normalization: $C=\epsilon$ frame

$S^3$ の構造定数は本稿全体で

$$
C^i{}_{jk}=\epsilon^{i}{}_{jk}\qquad(|C|=1;\ C=\epsilon\ \text{frame})
$$

に固定する。

**スケール変換と不変量.** 構造定数の全体スケール $C\to sC$（ $s>0$）の下で $n\to sn$ となるが、

$$
\mathrm{rank}(n),\qquad \dim\ker n,\qquad d_{\rm zero}=3-\mathrm{rank}(n)
$$

は不変である。一方、曲率と固有値の**数値**はスケール依存である。代表式:

$$
R\,q^2 \;\mapsto\; s^2\,(R\,q^2)\qquad(\text{等方点}),
$$

$$
\lambda_{S^3}(l)=(l-1)(l+3)\;\mapsto\; s^2\,(l-1)(l+3)\qquad(\text{TT sector; }l\ge2).
$$

本稿で $R\,q^2=6$、TT 固有値 $(l-1)(l+3)$ と書く箇所はすべて $C=\epsilon$（ $s=1$）frame の値である。`script/scripts/paper07/s3_normalization_factor_witness.py` で rank 不変と曲率・固有値の $s^2$ スケーリングを独立に確認できる。

**paper06 との対応**: paper06 [7] の $S^3$ 数値（例: $\chi=4$, $C_{\rm topology}=-36$, $V_3=2\pi^2$ ）は異なる構造定数正規化に基づく。両者は構造定数のスケール変換で関係づき、二次の曲率量には factor 16 の明示変換が入る。本稿では paper06 の $S^3$ 値を明示変換なしに引用しない。factor 変換は本 appendix に閉じ込める。

### A.4 Harmonic conventions

- $T^3$: Fourier モード $e^{ik\cdot x}$, $k\in 2\pi\mathbb Z^3/L$。
- $S^3$: scalar/vector/tensor harmonics は [28,29] の規約に従い、TT tensor harmonics の固有値は $C=\epsilon$ frame で $(l-1)(l+3)$（ $l\ge2$ ）と表す。文献 [28,29] の半径・正規化との対応は Appendix A.4 の規約に従う。
- $Nil^3$/ $Sol^3$: それぞれ toroidal ∪ Heisenberg sector、boost sector の取り方を [Appendix D](DPPUv7-paper07_appD.md) にまとめる。

### A.5 Diagnostic exclusions

本稿は Pontryagin 型 diagnostic（ $P_{\rm form}$, $P_{\rm int}$ 等 [5–7]）を導入しない。NY boundary 3-form の最小定義は [Appendix G](DPPUv7-paper07_appG.md) にのみ現れる（precision 目的）。

---

[← Section 9](DPPUv7-paper07_sec09.md) | [Appendix B →](DPPUv7-paper07_appB.md)
# DPPUv7-paper07 — Section 6: Harmonic One-Form Realization of the Zero Kernel

[← Section 5](DPPUv7-paper07_sec05.md) | [Section 7 →](DPPUv7-paper07_sec07.md)

---

## 6. Harmonic one-form realization of the zero kernel

本節では、Theorem 1 の scalar kernel が幾何学的対象としてどう実現されるかを示す。主張は left-invariant（代表的）スコープ内である。

### 6.1 Left-invariant harmonic one-forms and $\ker n$

> **Lemma 1（harmonic 1-form $\Leftrightarrow\ker n$）.**
> Bianchi class-A の left-invariant coframe $\{\sigma^i\}$ 上で、left-invariant harmonic 1-forms の空間は $\ker n$ と一致する。

**Proof.** left-invariant 1-form $\xi=\xi_i\,\sigma^i$（ $\xi_i$ 定数）に対し、構造方程式から

$$
d\xi = -\tfrac12\,\xi_i\,C^i{}_{jk}\,\sigma^j\wedge\sigma^k
$$

であり、class A では $C^{i}{}\_{jk}=\epsilon_{jkl} n^{li}$ なので、

$$
d\xi=0
\quad\Longleftrightarrow\quad
n\,\xi=0
\quad\Longleftrightarrow\quad
\xi\in\ker n.
$$

さらに left-invariant metric に対して co-closed 条件 $d{\star}\xi=0$ は unimodularity（class A）により自動的に満たされる。よって harmonic 1-forms $\cong\ker n$、次元は $d_{\rm zero}$。compact 商上では de Rham $H^1$ の代表元を与え、 $b_1$ との一致（[Section 4.5](DPPUv7-paper07_sec04.md)）の幾何学的内容を構成する [41,42]。外微分作用素 $d_1$（left-invariant 1-forms 上）と $n$ の nullspace 一致は `script/scripts/paper07/harmonic_kernel_correspondence.py` で確認可能である。

### 6.2 Map to the Hessian kernel

scalar-axis 摂動 $\beta_a$ から left-invariant 1-form 方向への写像を

$$
\Phi:\ \beta_a \longmapsto \beta_a\,\sigma^a
$$

（和なし; Milnor frame）と定義する。

> **Proposition 1（Hessian kernel $\leftrightarrow$ harmonic basis）.**
> reduced 二次作用の scalar-axis Hessian の kernel は $\Phi(\ker n)$ と一対一に対応する。 $\Phi$ は軸名ではなく subspace（次元 $d_{\rm zero}$）を保存する。

**Proof sketch.** $\beta\in\ker n$ では Lemma 1 により $\Phi(\beta)$ は harmonic。§4.3 の activation rank により $\beta$ 方向は secondary constraint を活性化せず、scalar-axis Hessian に kinetic pairing も発生しない。逆に Hessian kernel の各方向は $\ker n$ の基底 $\xi^{(\alpha)}$ に沿った frame 引き伸ばしとして実現される（[Appendix E](DPPUv7-paper07_appE.md)）。

### 6.3 Scope

この実現は left-invariant（homogeneous 代表）レベルの対応であり、非均質モードまで含む完全な Hodge 理論的主張ではない。非均質 sector で追加の zero が生じないことは [Section 5](DPPUv7-paper07_sec05.md) の非混入基準が独立に保証する。

---

[← Section 5](DPPUv7-paper07_sec05.md) | [Section 7 →](DPPUv7-paper07_sec07.md)
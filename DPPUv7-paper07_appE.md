# DPPUv7-paper07 — Appendix E: Harmonic Form and Hessian Kernel Realization

[← Appendix D](DPPUv7-paper07_appD.md) | [Appendix F →](DPPUv7-paper07_appF.md)

---

## Appendix E. Harmonic form and Hessian kernel realization

本 appendix では [Section 6](DPPUv7-paper07_sec06.md) の対応の検証表をまとめる。

### E.1 Lemma 1 and Proposition 1

[Section 6](DPPUv7-paper07_sec06.md) の Lemma 1（ $\mathrm{nullspace}(d_1)=\ker n$）と Proposition 1（写像 $\Phi:\beta_a\mapsto\beta_a\sigma^a$ が Hessian kernel と harmonic basis を対応）について、各 topology で以下を確認した:

1. **closed 条件**: left-invariant 1-form $\xi=\xi_i\sigma^i$ に対し $d\xi=0\Leftrightarrow n\xi=0$。外微分行列 $d_1$ と $n$ の一致 $\mathrm{nullspace}(d_1)=\mathrm{nullspace}(n)$。
2. **co-closed 条件**: class A（unimodular）の下で left-invariant 1-form は自動的に co-closed（ $d\star\xi=0$ ）。
3. **Hessian kernel 次元**: scalar-axis Hessian の kernel 次元 $=\dim\ker n$（軸名ではなく subspace 判定）。
4. **$\Phi$ の対応**: $\ker n$ の基底と Hessian kernel 方向の一対一（基底変換まで）。
5. **次元一致**: $\dim\{\text{harmonic 1-forms}\}=d_{\rm zero}=b_1$（ $T^3$:3, $Nil^3$:2, $Sol^3$:1, $S^3$:0）。

再現: `script/scripts/paper07/harmonic_kernel_correspondence.py`。

### E.2 Per-topology kernel directions（Milnor frame）

| topology | $n$ | $\ker n$ の基底 | harmonic 1-forms |
|---|---|---|---|
| $T^3$ | $\mathrm{diag}(0,0,0)$ | $e_1,e_2,e_3$ | $\sigma^1,\sigma^2,\sigma^3$ |
| $Nil^3$ | $\mathrm{diag}(1,0,0)$ | $e_2,e_3$ | $\sigma^2,\sigma^3$（fiber 方向 $\sigma^1$ は非調和） |
| $Sol^3$ | $\mathrm{diag}(1,-1,0)$ | $e_3$ | $\sigma^3$（単一特異軸） |
| $S^3$ | $\mathrm{diag}(1,1,1)$ | $\{0\}$ | なし |

軸名は Milnor frame 規約に依存し、frame 不変なのは次元のみである（[Section 4.4](DPPUv7-paper07_sec04.md)）。

### E.3 Scope note

この対応は left-invariant 代表レベルの realization である。compact 商上での de Rham $H^1$ 代表性は標準理論 [18] に依拠し、非均質 sector で追加の zero が生じないことは [Section 5](DPPUv7-paper07_sec05.md) の非混入基準が独立に保証する。

---

[← Appendix D](DPPUv7-paper07_appD.md) | [Appendix F →](DPPUv7-paper07_appF.md)
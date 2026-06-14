# DPPUv7-paper07 — Appendix G: Torsion Auxiliary Exactness and Field-Promotion Verification

[← Appendix F](DPPUv7-paper07_appF.md) | [Section 1 に戻る](DPPUv7-paper07_sec01.md)

---

## Appendix G. Torsion auxiliary exactness and field-promotion verification

本 appendix では [Section 8](DPPUv7-paper07_sec08.md) の検証詳細をまとめる。スコープは NY-MAIN route・AX/VT/MX torsion ansatz scalars である。

### G.1 Polynomial-degree verification（exact degree-2）

reduced torsion Lagrangian を [Section 8.2](DPPUv7-paper07_sec08.md) の

$$
\mathcal L_{\rm aux}= a_\eta\,\eta^2 + a_V\,V^2 + b_V\,\dot{\mathrm{vol}}\,V + \text{（boundary）}
$$

に整理し、 $(\eta,V)$ の多項式として次数 $D_{\max}=4$ まで symbolic に展開した。次数 3 以上の係数がすべて恒等的に 0（残差ゼロ）であることから、 $(\eta,V)$ について exact に高々二次である。

Euler–Lagrange 方程式は代数的に $\eta^\ast=0$、 $V^\ast=0$（ $V$ の線形項は $2\,d_t(\mathrm{vol}\cdot V)$ として境界に分離）を与え、打ち切りなしの **exact auxiliary shell** である。auxiliary Hessian（torsion mass datum）は負でなく、shell の安定性に整合する。

### G.2 Field-promotion derivative-coupling classification

$\eta(t)\to\eta(x)$, $V(t)\to V(x)$ の昇格に対し、EC+NY 作用に生じる微分結合を sector（EC/NY）× topology（4）× mode（AX/VT/MX）× 微分型（time/spatial/mixed）× 次数で分類した（計 128 cell）。結果:

| 微分型 | 二次（kinetic/gradient） | 線形 |
|---|---|---|
| time（ $\dot\eta,\dot V$ ） | **不在（全 cell）** | boundary 型 current のみ: EC $2\,\mathrm{vol}\cdot\dot V$ / NY $6\,\mathrm{vol}^{2/3}\,\dot\eta$ |
| spatial（ $\partial_i\eta,\partial_i V$ ） | **不在（全 cell）** | **線形項すら不在（全 cell）** |
| mixed | 不在 | 不在 |

二次微分結合は一切現れない。

**構造的根拠**: $\Gamma=\Gamma_{\rm LC}+K(\eta,V)$、 $R[\Gamma]=d\Gamma+\Gamma\wedge\Gamma$。 $K$ が $\eta,V$ に代数的なら $dK$ は field derivative に一次であり、 $\Delta R$ の scalar contraction では二次微分結合が構造的に消える（§8.3）。Riemann 成分レベル（24 成分）では個別に微分項が発火するが、scalar 縮約で相殺することを成分和の symbolic 計算で確認した。

**cross-check**: 同一の分類手続きを自由スカラー toy Lagrangian に適用し、kinetic 項を正しく検出することを確認した。再現: `script/scripts/paper07/torsion_auxiliary_field_promotion_witness.py`（full 24 成分 torsion 場理論への主張は含まない）。

### G.3 Minimal definition of the NY boundary 3-form（precision）

NY 4-form は本稿の規約で

$$
dQ_{\rm NY}=e_I\wedge e_J\wedge R^{IJ}-T^I\wedge T_I
$$

を満たす境界 3-form $Q_{\rm NY}$ を持つ（この恒等式は本稿の規約下で残差 0 を確認済み）。G.2 の NY 側時間線形 current $6\,\mathrm{vol}^{2/3}\dot\eta$ は $Q_{\rm NY}$ primitive の時間微分に一致し、boundary 性の整合をなす。 $Q_{\rm NY}$ はこの最小定義の範囲でのみ用い、本文の object としては導入しない（boundary charge の物理は本稿の範囲外）。

### G.4 Corollary table（finite fate invariance）

auxiliary shell 上で $L_{\rm aux}|_{\rm shell}=0$ であるため、[Appendix F](DPPUv7-paper07_appF.md) の modulus/wall 判定は EH/AX/VT/MX のすべてで不変である（系論的帰結; 独立力学検証ではない）。

### G.5 Scope and limitations (torsion sector)

- 本判定は AX/VT/MX ansatz scalars の場昇格に対するものであり、full 24 成分 torsion 場理論の $(\partial T)^2$ 不在定理ではない。
- 線形 boundary current を kinetic / propagating mode と呼ばない。
- exact degree-2 / exact shell は NY-MAIN route 限定（他 route へ自動拡張しない）。
- torsional-charge theorem / quantized charge / localization を主張しない。
- homogeneous auxiliary persistence（G.1）と field-promotion 検証（G.2）は別の層であり、後者は前者を再証明しない（拡張である）。

---

[← Appendix F](DPPUv7-paper07_appF.md) | [Section 1 に戻る](DPPUv7-paper07_sec01.md)
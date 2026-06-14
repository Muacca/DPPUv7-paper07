# DPPUv7-paper07 — Section 8: Torsion/Nieh–Yan Sector Robustness

[← Section 7](DPPUv7-paper07_sec07.md) | [Section 9 →](DPPUv7-paper07_sec09.md)

---

## 8. Torsion/Nieh–Yan sector robustness

作用には NY 項と捩れが含まれるにもかかわらず、ここまでの定理と頑健性検証は torsion に言及しなかった。本節では torsion sector が前述の結果に与える影響を検討する。主張は次の一点に限る: **torsion sector は、本稿の定理階層（二次 / $k>0$ スペクトル完成 / 有限変形）のいずれにも影響を与えない。** torsion 自体の物理（boundary charge・localization など）は本稿の範囲外である。

### 8.1 Known/new boundary

NY 項が正準解析において topological であり、古典力学を変えないことは既知である [33–35]。本節の内容はその一般論の再主張ではなく、**reduced 系における exact な実証**である: (i) reduced 作用が torsion ansatz scalars について厳密に高々二次であること、(ii) 場への昇格でも検査済み判定で微分結合が生じないこと、(iii) 有限変形の分類が不変であること。

### 8.2 (i) Exact auxiliary shell

NY-MAIN route の reduced torsion Lagrangian を auxiliary shell として

$$
\mathcal L_{\rm aux}
= a_\eta\,\eta^2 + a_V\,V^2 + b_V\,\dot{\mathrm{vol}}\,V + \text{（boundary total-derivative）}
$$

の形に整理できる（係数 $a_\eta,a_V,b_V$ は topology/mode 依存だが $\eta,V$ に関して代数的）。多項式次数検証（次数 4 まで symbolic、残差ゼロ）により、 $(\eta,V)$ について **exact に高々二次**である。

**代数的証明（shell 条件）.** Euler–Lagrange 方程式は

$$
\frac{\partial\mathcal L_{\rm aux}}{\partial\eta}=2a_\eta\,\eta=0,\qquad
\frac{\partial\mathcal L_{\rm aux}}{\partial V}=2a_V\,V + b_V\,\dot{\mathrm{vol}}=0
$$

を与える。 $a_\eta\ne0$ なら $\eta^*=0$。 $V$ の線形項 $b_V\,\dot{\mathrm{vol}}\,V$ は

$$
2\,\frac{d}{dt}\!\left(\mathrm{vol}\cdot V\right)
$$

として境界 current に分離でき、bulk 方程式は $V^*=0$ を与える。よって

$$
\eta^\ast=0,\qquad V^\ast=0
$$

は近似でない **exact auxiliary shell** である（[Appendix G](DPPUv7-paper07_appG.md)）。

### 8.3 (ii) Field promotion does not generate kinetic terms

接続を $\Gamma=\Gamma_{\rm LC}+K(\eta,V)$ と分解し、

$$
R[\Gamma]=d\Gamma+\Gamma\wedge\Gamma,\qquad
K=K(\eta,V)\ \text{（AX/VT/MX ansatz scalars に代数的）}
$$

と置く。 $K$ が $\eta,V$ に代数的なら $dK$ は field derivative に一次であり、 $\Delta R$ の scalar contraction における二次微分結合 $(\dot\eta^2,\dot V^2,(\nabla\eta)^2,(\nabla V)^2)$ は構造的に消える（native-action gate の主張; full torsion field theory ではない）。

scalars を場 $\eta(x), V(x)$ に昇格させ、上記を四 topology × AX/VT/MX の全 cell で検査した。結果:

- 二次微分結合は**一切現れない**。空間 gradient は線形項すら不在である。
- 出現するのは boundary 型の時間線形 current のみ: EC 側 $2\,\mathrm{vol}\cdot\dot V$、NY 側 $6\,\mathrm{vol}^{2/3}\,\dot\eta$。
- Riemann 成分レベルでは個別に微分項が発火するが、scalar 縮約で相殺することを成分計算で確認した。
- 分類手続きは自由スカラー toy Lagrangian の cross-check で kinetic 項を正しく検出することを確認済み（`script/scripts/paper07/torsion_auxiliary_field_promotion_witness.py`; full 24 成分 torsion 場理論への主張は含まない）。

### 8.4 (iii) Finite fate is unchanged

[Section 7](DPPUv7-paper07_sec07.md) の modulus/wall split は EH/AX/VT/MX のすべての torsion mode で不変である（auxiliary shell 上で $L_{\rm aux}|_{\rm shell}=0$ となることの系）。

### 8.5 Conclusion and scope

> NY-MAIN route と AX/VT/MX torsion ansatz scalars の範囲では、reduced torsion sector は exact auxiliary に留まり、 $\eta(x), V(x)$ への場昇格でも検査済み native-action 判定において二次 kinetic/gradient 項を生成しない。したがって torsion/NY は、DBR rank-defect kernel（[Section 4](DPPUv7-paper07_sec04.md)）・ $k>0$ スペクトル完成による頑健性（[Section 5](DPPUv7-paper07_sec05.md)）・finite homogeneous modulus/wall split（[Section 7](DPPUv7-paper07_sec07.md)）のいずれにも影響を与えない。

**適用範囲**: 本節の対象は AX/VT/MX ansatz scalars であり、full 24 成分 torsion 場理論に対する $(\partial T)^2$ 不在定理ではない。exact degree-2 / exact shell の主張は NY-MAIN route に限る。時間線形の boundary current を kinetic / propagating mode と呼ばない。NY の bulk への関与の不在を、このスコープの外で主張しない。

**Table 5** Torsion robustness の階層表.

| 階層 | 検査 | 結果 | スコープ |
|---|---|---|---|
| quadratic | 多項式次数検証（ $\eta,V$ ） | exact degree-2; $\eta^\ast=V^\ast=0$ exact shell | NY-MAIN route |
| finite | modulus/wall split の不変性 | 不変（系論的帰結） | finite Palatini branch |
| field-promoted | 微分結合分類（4 topology × AX/VT/MX） | kinetic/gradient 不在; boundary current のみ | 検査済み native-action 判定 |

---

[← Section 7](DPPUv7-paper07_sec07.md) | [Section 9 →](DPPUv7-paper07_sec09.md)
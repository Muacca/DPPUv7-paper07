# DPPUv7-paper07 — Appendix F: Finite-Fate Branch Data

[← Appendix E](DPPUv7-paper07_appE.md) | [Appendix G →](DPPUv7-paper07_appG.md)

---

## Appendix F. Finite-fate branch data

本 appendix では [Section 7](DPPUv7-paper07_sec07.md) の有限変形 branch データをまとめる。スコープは finite Palatini branch / homogeneous minisuperspace である。

### F.1 Branch parametrization and tangent alignment

kernel 方向ごとに、有限変形パラメータ $\beta$（および必要な補助方向 $u$ ）で background-adapted static Palatini branch を再構成する。branch の接ベクトルが二次解析の kernel 方向と整合すること（tangent alignment）を各 case で確認した。

### F.2 Energy-density expansion versus stationarity

- **energy-density 展開**: 補助場の方程式（auxiliary shell）を解いた上での reduced energy density $\Delta(\beta)$ の評価。本 appendix の主要検証。
- **停留点判定**: 有限変形点が field equations の停留点であるかの判定。

両者は異なる検査であり、energy-density 展開の結果（ $\Delta$ の値）を停留点判定の定理（stationary branch の存在）と読み替えない。

### F.3 Branch classification summary

**Branch 方向（Milnor frame; tangent alignment 確認済み）.**

| label | 接ベクトル $\boldsymbol\beta$ | 分類 |
|---|---|---|
| $T^3$ $e_a$ | $e_1,e_2,e_3$ | exact modulus |
| $Nil^3$ $\beta_-$ | $(\tfrac12,-\tfrac12,0)$ | exact modulus |
| $Nil^3$ $\beta_+/u$ | $(-\tfrac12,-\tfrac12,1)$ | curvature wall |
| $Sol^3$ $\beta_x$ | $(\tfrac13,-\tfrac16,-\tfrac16)$ | curvature wall |

**$\beta_-$ の closed-form branch 式（exact modulus の根拠）.** background-adapted static Palatini branch では scale factors が

$$
q_1(s)=q_0\,e^{s/2},\qquad q_2(s)=q_0\,e^{-s/2},\qquad q_3(s)=q_0
$$

（ $s$ は沿 $\beta_-$ の有限変形パラメータ、 $q_0$ は参照 scale）と書ける。接続 fluctuation を field equations で解いた reduced EH density $\Delta(s)$ は $s$ の全次数で恒等的に 0 である（4 次まで symbolic 残差 0 と一致）。 $\beta_+/u$ では $\Delta$ の 1 次係数（tadpole）が非零、 $\beta_x$ では線形・二次係数のいずれかが非零となる。

| 対象 | 検証 | 結果 |
|---|---|---|
| $T^3$（3 方向） | $\Delta(\beta)$ 4 次まで symbolic | 恒等的に 0（exact flat moduli） |
| $Nil^3$ $\beta_-$ | 同上＋線形系 rank 検証（rank = aug rank = 18, 一意可解） | 全次数残差 0（exact modulus） |
| $Nil^3$ $\beta_+/u$ | 同上 | 曲率由来の非零項＋tadpole（curvature wall; stationary 有限 branch なし） |
| $Sol^3$ $\beta_x$ | 同上 | 非零（curvature wall; stationary branch absent） |
| $S^3$ | — | 対象なし（ $d_{\rm zero}=0$ ） |

`script/scripts/paper07/finite_fate_branch_witness.py` は上記方向と分類を、finite Palatini branch / homogeneous minisuperspace のスコープ内で独立に再現する。最終的な branch 構成に基づく結果である。

### F.4 Minisuperspace closure verification

lever-complete ansatz（ $\log$-scale 変数 $l_a$ ・共役運動量 $p_a$ ・lapse $N$ ）の DBR 閉包:

- Hamiltonian constraint を含む閉包が成立（lapse 消去は閉包後）。
- DeWitt metric の trace-free tangent norm は正（minisuperspace no-ghost）。
- homogeneous anisotropy の configuration 自由度 = 2。

### F.5 Scope and limitations (finite-fate sector)

- wall は finite Palatini branch 内で確立されている。full nonlinear lift/mass theorem と呼ばない。
- energy-density 展開による分類を Hamiltonian-constraint theorem / 停留点判定と同一視しない。
- finite moduli 数（ $T^3$:3, $Nil^3$:1, $Sol^3$:0）と $d_{\rm zero}$（3, 2, 1）を同一視しない。
- minisuperspace no-ghost を full field-theory no-ghost と読み替えない。

---

[← Appendix E](DPPUv7-paper07_appE.md) | [Appendix G →](DPPUv7-paper07_appG.md)
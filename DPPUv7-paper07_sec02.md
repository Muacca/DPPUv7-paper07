# DPPUv7-paper07 — Section 2: Lorentzian EC+NY Setup and Bianchi Class-A Backgrounds

[← Section 1 (Introduction / References)](DPPUv7-paper07_sec01.md) | [Section 3 →](DPPUv7-paper07_sec03.md)

---

## 2. Lorentzian EC+NY setup and Bianchi class-A backgrounds

### 2.1 First-order action ingredients

基本変数は tetrad 1-form $e^I$ と Lorentz 接続 1-form $\omega^{IJ}$ であり、両者を独立変数として扱う（first-order / Palatini 型）。作用は Einstein–Cartan 項 [8] と Nieh–Yan 項 [9] からなる:

$$
S[e,\omega]
=\frac{1}{2\kappa^2}\int \epsilon_{IJKL}\, e^I\wedge e^J\wedge R^{KL}[\omega]
+\frac{\theta_{\rm NY}}{\kappa^2}\int \left( T^I\wedge T_I - e_I\wedge e_J\wedge R^{IJ} \right),
$$

ここで $R^{IJ}[\omega]=d\omega^{IJ}+\omega^{I}{}_K\wedge\omega^{KJ}$、 $T^I=de^I+\omega^I{}_J\wedge e^J$ である。signature・Hodge/ $\epsilon$ 規約は [Appendix A](DPPUv7-paper07_appA.md) に固定する。

**NY 項の取り扱い**: 本稿では NY 項を exact 4-form として扱う route（以下 NY-MAIN route）を採用する。**本稿の NY-MAIN route と LZ-native AX/VT/MX torsion ansatz scalars の検査済み reduced/native-action gate では**、背景上の NY 項は exact boundary term としてのみ寄与する。reduced 系では exact auxiliary shell をなし、場 $\eta(x), V(x)$ への昇格でも検査済み native-action 判定において二次 kinetic/gradient 項は生じない（[Section 8](DPPUv7-paper07_sec08.md)）。NY 項の topological 性は正準解析で既知である [33–35]。**この bulk 無関与の主張は上記スコープ内に限定**し、full torsion field theory や一般の NY 作用に自動拡張しない。

### 2.2 Bianchi class-A backgrounds and the structure matrix

空間断面は 3 次元 unimodular Lie 群 $G$ の left-invariant 構造を持つ。left-invariant 1-forms $\sigma^i$（ $i=1,2,3$ ）は

$$
d\sigma^i=-\tfrac12\,C^i{}_{jk}\,\sigma^j\wedge\sigma^k
$$

を満たす。class A（unimodular）条件は $C^j{}_{ij}=0$ であり [15]、このとき構造定数は対称な**構造行列** $n^{ij}$ で表される:

$$
C^i{}_{jk}=\epsilon_{jkl}\,n^{li},
\qquad
n^{ij}=\tfrac12\,\epsilon^{ikl}\,C^{j}{}_{kl},
\qquad n^{ij}=n^{ji}.
$$

Milnor frame [16] により $n^{ij}=\mathrm{diag}(n_1,n_2,n_3)$ と対角化でき、その rank/signature が class-A 型を分類する。本稿で扱う四 topology を Table 1 にまとめる。compact な空間断面を持つ Bianchi モデルは class A に限られる [19,20]。

**Table 1** Four Bianchi class-A backgrounds.

| topology | Bianchi type | $n$ (Milnor frame) | $\mathrm{rank}(n)$ | $d_{\rm zero}=3-\mathrm{rank}(n)$ | $b_1$ | $R\,q^2$（等方点） |
|---|---|---|---:|---:|---:|---:|
| $T^3$ | I | $\mathrm{diag}(0,0,0)$ | 0 | 3 | 3 | $0$ |
| $Nil^3$ | II | $\mathrm{diag}(1,0,0)$ | 1 | 2 | 2 | $-1/2$ |
| $Sol^3$ | VI₀ | $\mathrm{diag}(1,-1,0)$ | 2 | 1 | 1 | $-2$ |
| $S^3$ | IX | $\mathrm{diag}(1,1,1)$ | 3 | 0 | 0 | $6$ |

ここで $b_1$ は対応する compact 商の first Betti number（標準的な値 [18]）、 $R\,q^2$ は等方点での無次元化された空間スカラー曲率である。 $S^3$ の構造定数正規化は $C=\epsilon$（ $|C|=1$ ）に固定する（[Appendix A](DPPUv7-paper07_appA.md); paper06 [7] の正規化との明示変換もそこに記す）。 $\mathrm{rank}(n)$ と $d_{\rm zero}$ は構造定数の全体スケール $C\to sC$ に対して不変だが、 $R\,q^2$ や TT 固有値などの**数値**は $s$ に依存する（式は Appendix A.3）。

Table 1 の $d_{\rm zero}$ は構造行列の rank-defect による代数的予測である。これが制約スペクトルの kernel として実現されることが本稿の主結果である（[Section 4](DPPUv7-paper07_sec04.md)）。abelianization との関係 $\dim(\mathfrak g/[\mathfrak g,\mathfrak g])=\dim\ker n$ も既知の代数的事実であり、幾何学的 cross-check としてのみ用いる。

### 2.3 Background ansatz and perturbation sectors

背景は homogeneous な対角 ansatz

$$
e^0=N(t)\,dt,\qquad e^a=q_a(t)\,\sigma^a \quad (a=1,2,3;\ \text{和なし})
$$

をとり、接続 $\omega^{IJ}$ は独立な homogeneous 成分を持つ（具体的な成分リストと reduced action は [Appendix B](DPPUv7-paper07_appB.md)）。摂動は背景 frame に関して

- **scalar（diagonal/axis）sector**: 対角 metric 摂動 $\beta_a$ と対応する接続摂動、
- **spin-1（transverse vector）sector**,
- **spin-2（TT）sector**: transverse-traceless 摂動 $h_{ab}^{\rm TT}$,

に分解する。homogeneous モード（ $k=0$ ）の解析が [Section 4](DPPUv7-paper07_sec04.md) の中心であり、非均質モードは topology ごとの harmonic 展開（ $T^3$: Fourier、 $S^3$: hyperspherical harmonics [28,29]、 $Nil^3$: toroidal ∪ Heisenberg sector、 $Sol^3$: boost sector）で扱う（[Section 5](DPPUv7-paper07_sec05.md), [Appendix D](DPPUv7-paper07_appD.md)）。

### 2.4 Torsion ansatz scalars

homogeneous 背景上の torsion は、本稿で用いる ansatz の範囲で次の scalar modes に整理される:

- **AX（axial）**: $T^{a}{}\_{bc}\supset \eta\,\epsilon_{abc}$ を与える axial scalar $\eta(t)$、
- **VT（vector/trace）**: trace part を与える scalar $V(t)$、
- **MX（mixed）**: 両者の混合 branch。

これらは paper05/06 [6,7] の branch 分類（EH/AX/VT/MX）と同一の系統である。[Section 8](DPPUv7-paper07_sec08.md) ではこれらを場 $\eta(x), V(x)$ に昇格させた際の微分結合を分類する。本稿の torsion に関する主張はすべて **AX/VT/MX ansatz scalars の範囲内**であり、full 24 成分 torsion 場理論への一般化は主張しない。

### 2.5 Notation summary

以降、 $q^2$ は等方 scale 因子の二乗（ $q_1=q_2=q_3=q$ の等方点）、 $\kappa^2=8\pi G$、 $\mathrm{vol}$ は空間 cell の comoving 体積を表す。制約解析の記号（primary/secondary constraints, first/second class, Dirac bracket）は標準的な定義 [10–12] に従う。

---

[← Section 1](DPPUv7-paper07_sec01.md) | [Section 3 →](DPPUv7-paper07_sec03.md)
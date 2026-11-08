# 6. 非斉次摂動の拘束還元

## 共通の拘束・時間依存・実場規約

各blockで六metric成分、lapse、三shiftを保持して作用を変分する。複素表示では共役sectorを最後に復元し、非対角momentumは $`p_{ij}=2\pi^{ij}`$ とする。z=(h,p)、Ω=−J、gauge行列G、埋込みEを使うと、共通の還元式は次となる。

式(R-01)。

$$
\begin{aligned}
CJC^\dagger&=0,\qquad(\dot C+CJH)\ker C=0,\qquad CE=GE=0,\\
E^\dagger\Omega E&=\Omega_r,\qquad H_r=E^\dagger HE-E^\dagger\Omega\dot E,\\
\dot y&=J_rH_ry,\qquad z=Ey.
\end{aligned}
$$

随伴は必要に応じ空間部分積分を含む。時間依存する拘束と埋込みの微分は省略しない。還元作用だけでなくmultiplierを復元して元の方程式へ戻す。付録F–Gに基底・閉包・gauge・非零minor・再構成を示す。

## 固定flat torus

N₀=1、g₀=q₀²Iは完全真空背景である。計数ではc=q₀=1、非零波数をz方向に表示し、平均二乗1の実基底を用いる。metricとlapseにcos、shiftにsinを使い、未還元二次作用を作る。

式(W-01)。

$$
\begin{aligned}v_{ij}&=\dot h_{ij}-k(\beta_i\delta_{j3}+\beta_j\delta_{i3}),\\
L^{[2]}_k&=\tfrac14\bigl[\mathrm{tr}(v^2)-(\mathrm{tr}v)^2\bigr]+\mathcal R^{[2]}(h)+\alpha k^2(h_{11}+h_{22}),\\
\mathcal R^{[2]}&=\delta^{ij}(\Gamma^{k(1)}{}_{i\ell}\Gamma^{\ell(1)}{}_{jk}-\Gamma^{k(1)}{}_{ij}\Gamma^{\ell(1)}{}_{k\ell}).\end{aligned}
$$

四primaryと四secondaryは第一級である。secondaryはp13=p23=p33=0とh11+h22=0に比例し、gauge商に二対が残る。波数の表示を変えるframe回転は格子の自己同型を仮定しない。k=0には非零kでの割算を適用せず、元の二次作用で六metric pairへ戻る。

式(W-02)。

$$
\begin{aligned}h_{11}&=A,\quad h_{22}=-A,\quad h_{12}=B,\quad h_{i3}=0,\\
L^{[2]}_{\rm TT}&=\tfrac12\bigl(\dot A^2+\dot B^2-k^2(A^2+B^2)\bigr),\\
S^{[2]}_{\rm TT}&=\frac{\mathcal V_{\Sigma}q_0^3}{8\kappa^2}\int d\tau\thinspace \left\langle\dot\gamma_{ij}\dot\gamma_{ij}-q_0^{-2}\partial_a\gamma_{ij}\partial_a\gamma_{ij}\right\rangle.\end{aligned}
$$

γは平均零のTT場である。一つの実Fourier成分の二偏極は正kineticで、静的torusではω²=|k|²/q₀²となる。完全な実±k blockの総数は4対である。この計数と第7節の厳密族への積分条件を区別する。

## Nil：central零と非零sector

central指数m=0でも、base波数(n₁,n₂)≠0なら非斉次である。固定LRS背景ではcovectorをframe内で(k,0)へ揃え、k²=n₁²+n₂²とする。格子を回転・同定し直す操作ではない。R=4u/3、ζ=16u²/[9(1+u⁴)²]により還元作用は次となる。

式(N-01)。

$$
L^{[2]}=\left\langle\frac R6(P_u^2+\zeta Q_u^2)
-\frac{3R}{2}\bigl[P_x^2+\zeta(Q_x^2-4PQ_x+2P^2)\bigr]\right\rangle.
$$

微分結合はcos/sinを混ぜる。P∝cos、Q∝sinの実blockとその相補blockに各2対、完全な実反対波数blockには4対ある。採用区間のζ≤4/9とk²≥1からkineticとpotentialは正だが、Hamiltonianは時間依存し、保存energyや長時間安定性を意味しない。

m≠0ではHermite–theta scalar関数に全tensor spin sを組み合わせる。LRS係数とladder則によりL=n+sign(m)sが保存され、lapse・shift・拘束・gaugeを含む有限blockが閉じる。Heisenberg表現自体は無限次元であり、固定商上の回転Killing場を仮定しない。

| L | 複素metric | 拘束rank | 複素対 |
| --- | --- | --- | --- |
| -2 | 1 | 0 | 1 |
| -1 | 2 | 1 | 1 |
| 0 | 4 | 3 | 1 |
| 1 | 5 | 4 | 1 |
| 2以上 | 6 | 4 | 2 |

表は固定m>0・格子コピーrの計数で、共役−mを含む実場では倍、全rではさらにm倍する。低端ではn<0のtensor・multiplierが存在せず、残る拘束は全行rankを持つ。一般L≥2の計数は閉包と非零central KIDがないことから従う。明示行列はL=−2,−1,0,1,2、独立4次元照合はL=2で与える。非LRSではn↔n+2の結合が残り、同じ有限gradingへは閉じない。

## Sol：Gowdy sectorと非零dual orbit

fiber非依存sectorでは、P=2x+p、Q=e⁻²ˣrを用いるとp,rはbase周期 $`L_{\rm Sol}=\log[(3+\sqrt5)/2]`$ を持つ。非零の $`k=2\pi n/L_{\rm Sol}`$ に対して、実cosまたはsin成分の還元作用は次となる。

式(SO-01)。

$$
L_k^{[2]}=\frac t2\left[\dot p^2+\dot r^2-k^2p^2-(k^2+4)r^2\right].
$$

これは二つの伝播場で、完全な実blockには4対ある。k²+4はareal時間と選んだ変数の係数であり、保存物理振動数ではない。一般の非零fiber sectorではdual latticeのmonodromy orbit全体を保持し、係数を実数直線上の関数として表す。

式(SO-02)。

$$
D_1=\partial_\xi,\qquad D_2=iK,\qquad D_3=iM,
\qquad K=k_0e^{-\xi},\quad M=l_0e^\xi.
$$

全metric・lapse・shift・拘束・gaugeは同じ無限orbitを保つ。固定Hermite基底では指数関数の乗算が無限個の係数を結ぶ。付録Gの微分右逆は、固定格子の全非零orbitで分母が零から離れることを利用し、二つの複素正準場対へ還元する。実場には共役orbitを合わせて4実正準場対がある。これは無限次元の場の種類数であり、有限blockの総数2ではない。

## S³ Taub：完全なWigner block

不変frame成分をWigner関数に展開し、固定jの全右添字を保持する。左重複度と実場条件を復元した結果は次の通りである。jはround sphereのTT harmonicラベルとは異なる。

| j | 実拘束rank | 実正準対 | KID欠損 |
| --- | --- | --- | --- |
| 0 | 3 | 3 | 軸U(1) |
| 1/2 | 16 | 8 | なし |
| 1 | 33 | 21 | 左SU(2) |
| 3/2 | 64 | 32 | なし |

j=0、1の値は対称背景の線形gauge商であり、二次KID条件が残る。Nil低端の成分不在とは異なる機構である。j=3/2では二場に各16実Wigner成分が属し、32対となる。gauge、時間依存埋込み、multiplierを用いた再構成は元の4次元真空方程式の全成分を満たす。表は計算した四ラベルだけを示し、全jの一様分類は主張しない。

## 有限時間例と独立なSol対照

付録Hの指定初期値では、Sol GowdyのBessel評価、Nil baseの一階・二階発展、S³ j=3/2およびNil m=1,r=0,L=2の元ADM発展との比較が整合する。これらは非自明な線形例であり、曲がった背景の非線形持ち上げを示すものではない。

Solの一orbit (1,0)では、同じ初期ADMデータから独立に導いた四次元Lorenz系を発展させた。1≤t≤2、内部領域|ξ−ξ₀|≤4で、還元変数のL²差は1023→1535点で3.26×10⁻⁸、区間幅対照で2.77×10⁻¹⁰、時間対照で3.94×10⁻⁹だった。実際の時間補間を微分した全10 Einstein成分のL²残差は20検査時刻で最大3.83×10⁻¹¹である。

一方、再構成した還元変数を4階還元系SG-05へ代入した正規化強残差（SH-05）は、同じ検査で最大1.54×10⁻³である。四つの軌道の保存時刻で評価すると9.6×10⁻⁴から4.5×10⁻³の範囲にあり、1023点から1535点へ細かくしても減らず、時間許容誤差を緩めると大きくなる。したがって、この形での収束は確立していない。二つの残差は異なるnormである。有限区間のDirichlet条件、内部の特性余裕、区間幅対照を付録Hに記す。この例を実数直線全体の誤差保証や一般安定性へ拡張しない。

内部profileの図を[付録H](DPPUv7-paper07_appH.md)に示す。

[論文の入口](README_ja.md) / [参考文献](REFERENCES.md)

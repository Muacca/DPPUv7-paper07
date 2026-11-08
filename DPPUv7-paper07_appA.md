# 付録A. 固定商と規約の辞書

T³の座標x,y,zは周期2πとし、σ=(dx,dy,dz)を用いる。Nilではσ=(dx,dy,dz−xdy)とし、次の三つの格子生成元を固定する。

式(QA-01)。

$$
(x,y,z)\sim(x+2\pi,y,z+2\pi y)\sim(x,y+2\pi,z)\sim(x,y,z+4\pi^2).
$$

各生成元はcoframeを保ち、基本領域の体積は(2π)⁴である。central Fourier因子はexp(imz/(2π))、mは整数である。双対場X=∂x、Y=∂y+x∂z、Z=∂zから−(X²+Y²+Z²)を作ると、m=0はbase torus、m≠0は周波数|m|/(2π)のoscillatorとcentral定数m²/(4π²)になる。twisted periodicityは符号付きmごとに|m|個の独立状態を与える。

Solではσ=(dx,eˣdy,e⁻ˣdz)、A=((2,1),(1,1))、λ=(3+√5)/2とする。Aの固有座標で面積1のfiber格子を取り、(x,y,z)と(x+log λ,λ⁻¹y,λz)を同定する。構造行列はn₂₃=n₃₂=1、参照体積はlog λである。fiber Fourier波数はdual latticeのA作用のorbitとしてまとめる。

S³はsimply connected SU(2)で、dσⁱ=2σʲ∧σᵏ（循環）を満たす単位round coframeを用いる。q²Σσⁱ²の断面曲率は1/q²、scalar曲率は6/q²、体積は2π²q³である。構造方程式を負号で定義する場合は、nの符号を反転して本稿の規約へ換算する。一般の軸変換は向きを含めた同じ行列をmetric・接続・torsionにも作用させる。

式(QA-02)。

$$
C^i{}_{jk}=\frac{a}{r}\epsilon_{ijk}\quad\Longrightarrow\quad q=\frac{2r}{|a|},\quad R=\frac{3a^2}{2r^2},\quad {\rm Vol}=\frac{16\pi^2r^3}{|a|^3}.
$$

式(QA-02)は物理正規直交frameでの構造係数の大きさをa/rとする換算である。a=1ならq=2r、a=2ならq=r、a=4ならq=r/2となる。曲率だけを換算して体積や波数を据え置かない。本稿では幾何学的torsion Tを直接使用する。実装間の換算は計算補足の比較ガイドに示す。

式(QA-03)。

$$
\mathcal N=\frac14[abcd]T^e{}_{ab}T_{ecd}-\frac12[abcd]R_{abcd},\qquad[0123]=+1.
$$

この交代記号は、Lorentz計量で添字を上げたepsilonテンソルではない。完全NYの両項を同一の接続から作り、式(A-01)の外微分と照合する。

## 摂動の規格化と記号

以下の非斉次計算では正の全体作用係数を除き、固定Haar測度で積分する。対称metricの順序は(11,22,33,12,13,23)、共役運動量は対角で $`\delta\pi^{ii}`$ 、非対角で $`2\delta\pi^{ij}`$ とする。この規約でJは標準正準行列となる。daggerは複素共役と、空間積分に関する形式随伴を含む。

式(QA-04)。

$$
z=(h,p),\qquad J=\begin{pmatrix}0&I\\
-I&0\end{pmatrix},\qquad
\Omega=-J,\qquad [X_i,X_j]=C^k{}_{ij}X_k.
$$

Nilのmはcentral Fourier整数、rは格子コピー、nはoscillator次数、Lはtensor gradeである。S³のj,m,nはWignerラベルであり、背景のTaub parameter m=0とは独立である。Solでは元のbase座標をx、unfolding後の実数座標をξとする。以下でPは文脈によりNil周期2πまたはGowdy場を表す。混同を避けるため周期はP₀とも記す。NG-02のLRS背景値は $`\bar P`$ 、 $`\zeta=e^{2\bar P}`$ と記す。Solの周期は $`L_{\rm Sol}=\log\lambda`$ と記し、作用の $`L_0`$ 、 $`L^{[2]}`$ やNilのgrade Lと区別する。

[論文の入口](README_ja.md) / [参考文献](REFERENCES.md)

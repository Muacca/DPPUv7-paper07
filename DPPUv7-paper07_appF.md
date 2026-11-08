# 付録F. 真空背景・固定商・調和基底

## 背景の導出と商への降下

BF-01のNilはTanimotoの斉次真空解のb=1、Kasner指数(2/3,2/3,−1/3)を取り、文献時間をu³へ換えたものである。Solはtwisted Gowdy表示のP=2x,Q=0,R=tであり、S³はbiaxial Taubのℓ=1,m=0をσ=−ω/2へ換算した。全4次元Ricci tensor、Hamiltonian・momentum拘束、同じlapseによるmetricと運動量の発展を直接代入して零を得る。時間変換は固定格子を変えない。 背景の一次資料は[Tanimoto2004](REFERENCES.md#tanimoto2004)、[Rendall2012](REFERENCES.md#rendall2012)、[BeyerHennig2012](REFERENCES.md#beyerhennig2012)である。

式(BF-02)。

$$
\begin{aligned}
\pi_{\rm Nil}^{ij}&=\mathrm{diag}(-4/(9A),-4/(9A),-u^2),\\
\pi_{\rm Sol}^{ij}&=\mathrm{diag}(-1/F,-t-1/(4t),-t-1/(4t)),\\
\pi^{ij}&=\sqrt{\det G}\thinspace (K^{ij}-G^{ij}K),\qquad K_{ij}=\dot G_{ij}/(2N).
\end{aligned}
$$

Nilの三deck生成元はQA-01のcoframeをそのまま保つ。Solの貼り合わせも(eˣdy,e⁻ˣdz)を保ち、S³のcoframeは大域的である。したがって時間だけに依存するGとNは全て固定商に降下する。BF-01の閉区間で正定値性は表示から直接従う。

## Nilのcentral非零基底とtensor grade

式(NF-01)。

$$
\begin{aligned}
\mu&=m/P_0,\quad m\in\mathbb Z\setminus\lbrace 0\rbrace ,\quad 0\le r<|m|,\quad P_0=2\pi,\\
\Phi_{m,r,n}&=\sqrt{P_0}\thinspace e^{i\mu z}\sum_{k\in\mathbb Z}
 e^{i(r+km)y}f_n\left(x+\frac{r+km}{\mu}\right),\\
f_n(v)&=\frac{(|\mu|/\pi)^{1/4}}{\sqrt{2^nn!}}
 H_n(\sqrt{|\mu|}v)e^{-|\mu|v^2/2}.
\end{aligned}
$$

測度をdx dy dz/P₀⁴に規格化する。格子のxシフトは和の添字をずらし、y,z周期は整数位相を与えるためNF-01は商に降下する。異なるrとmは直交し、x区間の和を実数直線へ展開するとfₙの正規直交性からΦの正規直交性が従う。

式(NF-02)。

$$
\begin{aligned}
(X-iY)\Phi_n&=\sqrt{2\mu n}\thinspace \Phi_{n-1},&
(X+iY)\Phi_n&=-\sqrt{2\mu(n+1)}\thinspace \Phi_{n+1},&Z\Phi_n&=i\mu\Phi_n\quad(m>0),\\
\theta^\pm&=(\sigma^1\pm i\sigma^2)/\sqrt2,&\theta^0&=\sigma^3,& L&=n+\mathrm{sign}(m)s.
\end{aligned}
$$

対称tensor基底はθ⁺θ⁺、(θ⁺θ⁰+θ⁰θ⁺)/√2、(θ⁺θ⁻+θ⁻θ⁺)/√2、θ⁰θ⁰、負spinの二つで、spinは2,1,0,0,−1,−2である。m>0ではn=L−s≥0の成分だけを残す。LRS背景はspin零であり、微分のladder移動とtensor spin移動がLを保存する。metric、momentum、lapse、shift、Cartanの各tensor型を同じ規則で構成すると有限閉包が得られる。m<0は複素共役であり、rは−r modulo |m|、spinは−sへ移る。

このspinは基底を整理する局所の回転gradingであり、固定商の新しい回転KIDではない。非LRS係数はnとn±2を結合し、この有限grade閉包を壊す。非LRSでの解析そのものが不可能という意味ではない。m=0では通常のbase Fourier基底を用いる。

## 固定Sol格子とunfolding

式(SF-01)。

$$
\begin{aligned}
\varphi&=(1+\sqrt5)/2,\quad\lambda=\varphi^2,\quad L_{\rm Sol}=\log\lambda,\\
S&=\frac1{\sqrt{1+\varphi^2}}\begin{pmatrix}1&-\varphi\\
\varphi&1\end{pmatrix},\quad
D=\mathrm{diag}(\lambda^{-1},\lambda),\quad SA=DS,\\
q&=2\pi S(n_1,n_2)^T=(k_0,l_0)^T,\qquad q_j=(D^T)^jq,\\
k_0l_0&=\frac{4\pi^2}{\sqrt5}(n_1^2-n_1n_2-n_2^2)\ne0\quad(n\ne0).
\end{aligned}
$$

fiber格子はSℤ²で面積1である。SF-01の整数二次形式は非零整数対で零にならず、全ての非零dual orbitに両fiber成分が存在する。軸上の連続covectorをこの固定商のmodeに含めない。

式(SF-02)。

$$
\begin{aligned}
U_qf&=\sqrt{L_{\rm Sol}}\sum_{j\in\mathbb Z}f(x+jL_{\rm Sol})
 e^{i(\lambda^{-j}k_0y+\lambda^jl_0z)},\\
D_1&=\partial_\xi,\quad D_2=iK,\quad D_3=iM,
\qquad K=k_0e^{-\xi},\quad M=l_0e^\xi,\\
\mathcal D&=\lbrace f:e^{a|\xi|}\partial_\xi^bf\in L^2(\mathbb R) \text{for all integers }a,b\ge0\rbrace .
\end{aligned}
$$

正規化測度 $`dx\,dy\,dz/L_{\rm Sol}`$ ではUはL²(R)からorbit空間への等長写像となる。deck変換の添字移動とfiber直交性から直接確認できる。D上で全係数・微分・形式随伴を扱い、境界項は零である。compact supportを持つ滑らかな関数とHermite関数はDに含まれる。この形式的共通定義域の構成は閉作用素拡張やwellposedness定理ではない。

式(SF-03)。

$$
[D_1,D_2]=-D_2,\quad[D_1,D_3]=D_3,\quad
e^\xi\psi_0= e^{1/4}\sum_{n=0}^{\infty}\frac{\psi_n}{2^{n/2}\sqrt{n!}}.
$$

SF-03の全係数が非零なので有限Hermite切断は厳密には閉じない。無限orbit上では全tensor成分、拘束、gauge、発展が閉じる。共役orbitを合わせて実場とする。一方fiber零ではx周期のbase Fourierが適切で、両sectorを区別する。

## S³のWigner blockと実場

式(TF-01)。

$$
\begin{aligned}
Y^j_{mn}&=\sqrt{2j+1}\thinspace D^j_{mn},\quad d=2j+1,\quad
\overline{Y^j_{mn}}=(-1)^{m-n}Y^j_{-m,-n},\\
D_a&=2iJ_a,\quad[D_a,D_b]=-2\epsilon_{abc}D_c,\quad
\sum_aD_a^2=-4j(j+1)I.
\end{aligned}
$$

Haar測度を体積で割って正規化する。背景coframeに関する六metric係数をそれぞれ全右添字で展開する。左添字は重複度であり、固定jの完全な実blockはd²実scalar成分を持つ。tensor係数もTF-01と同じ共役条件を満たす。j=0,1/2,1,3/2だけを明示分類し、j=3/2を完全再構成例とする。

[論文の入口](README_ja.md) / [参考文献](REFERENCES.md)

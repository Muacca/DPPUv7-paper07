# 2. 作用・幾何・境界

signatureは $`(-,+,+,+)`$ 、向きは正の $`e^{0123}`$ 、 $`\epsilon_{0123}=+1`$ とする。coframeは非退化で空間sliceはspacelike、独立接続はmetric-compatibleなLorentz接続である。NY結合θは定数、κ²は正とし、物質・spin source・宇宙項・曲率二乗項を含めない。

式(A-01)。

$$
\begin{aligned}T^a&=de^a+\omega^a{}_b\wedge e^b,\qquad R^a{}_b=d\omega^a{}_b+\omega^a{}_c\wedge\omega^c{}_b,\\
Q&=e^a\wedge T_a,\qquad dQ=T^a\wedge T_a-e^a\wedge e^b\wedge R_{ab},\\
S&=\frac1{4\kappa^2}\int\epsilon_{abcd}e^a\wedge e^b\wedge R^{cd}+\theta\int dQ+S_{\partial}.\end{aligned}
$$

全ての反対称添字を和するとepsilon外積は $`2R_{\rm EC}e^{0123}`$ に等しい。したがって最初の項は $`R_{\rm EC}/(2\kappa^2)`$ となる。この規格化を正準運動量とkinetic係数にも一貫して用いる。

式(G-01)。

$$
d\sigma^i=\frac12 n_{ij}\epsilon_{jkl}\sigma^k\wedge\sigma^l,\quad n^T=n,\quad \mathcal V_{\Sigma}=\int_{\Sigma}\sigma^{123},\quad g=g_{ij}(t)\sigma^i\sigma^j.
$$

nは固定したclass Aの構造行列で、gは任意の正定値対称行列である。一定SO(3)回転はnのrankを保つ。本稿の固定商と参照規格化は次表のとおりで、球面のqは物理半径である。商の具体的同定は付録Aに示す。

| 幾何 | nの固有値 | 参照体積 | R q² |
| --- | --- | --- | --- |
| T³ | (0,0,0) | $`(2\pi)^3`$ | 0 |
| Nil | (0,0,-1) | $`(2\pi)^4`$ | $`-1/2`$ |
| Sol | (0,1,-1) | $`\log((3+\sqrt5)/2)`$ | $`-2`$ |
| S³ | (2,2,2) | $`2\pi^2`$ | 6 |

## 全接続方程式

式(A-02)。

$$
\epsilon_{abcd}e^c\wedge T^d=0\quad\Longrightarrow\quad T^a=0,\qquad \omega=\omega_{\rm LC}[e].
$$

定数完全NYはこのbulk方程式に寄与しない。epsilonを逆変換するとeᵃ∧Tᵇ−eᵇ∧Tᵃ=0となる。t=ΣᵃιₐTᵃと置いて内部積を取ればTᵇ=−eᵇ∧t、さらにt=−3tを得るためt=T=0である。AX/VT/MXへの制限を使わない証明である。付録Bには全24成分による代数的対照を与える。

K=ω−ωLCとすると、EC scalarはLC scalar、Kの代数的二次形式、発散に分かれる。発散を端点項へ移すとKは可逆な補助場となる。ωLCがcoframe微分に依存するため、単なるpoint変換として正準同値を仮定せず、Euler chain ruleと端点1形式を使う。

## 境界代表

時間両端のspacelike sliceでは、外向き法線の符号とε_n=−1を含むLC GHY項を加える。主線のbulk代表S₀は、さらにcontortion発散とNY currentを端点から差し引いたもので、固定誘導計量に対する変分に用いる。currentを保持した代表は、端点currentの変分を固定する対応データで比較する。境界sourceを導入しない。

式(A-03)。

$$
\begin{aligned}S_{\rm GHY}&=\kappa^{-2}\sum_{i,f}\epsilon_n\int_{\Sigma}\sqrt h\thinspace K_{\partial},\qquad \epsilon_n=-1,\\
S_{\rm keep}&=S_0+[B_{\rm end}]_i^f,\quad \Theta_{\rm keep}=\Theta_0+\delta B_{\rm end},\quad \Omega_{\rm keep}=\Omega_0.\end{aligned}
$$

共通scale・二振幅族では、幾何学的torsionをTⁱ=(V/3)e⁰ⁱ+(2η/q)eʲᵏとする。D=N⁻¹d/dt、H=Dq/qとおくと、以下の発散と端点関数になる。HVだけを全微分として削除してはならない。

式(A-04)。

$$
\begin{aligned}R_{\rm EC}-R_{\rm LC}&=-6\eta^2/q^2+2V^2/3-2DV-6HV,\\
B_{\rm end}&=-\mathcal V_{\Sigma}q^3V/\kappa^2+6\theta\mathcal V_{\Sigma}q^2\eta,\\
P_i&=p_i+\partial_iB_{\rm end}.\end{aligned}
$$

この代表で端点のq、η、Vを固定すればδB_end=0となる。lapseやscale速度を固定する必要はない。滑らかな完全接続shellでは端点でもη=V=0である。発散の除去はkinetic項の存在を変更せず、時間線形のcurrentは伝播モードではない。

[論文の入口](README_ja.md) / [参考文献](REFERENCES.md)

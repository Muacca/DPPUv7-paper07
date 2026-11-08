# 付録D. Fourier作用と積分拘束の導出

式(WD-01)。

$$
\Gamma^{i(1)}{}_{jk}=\tfrac12(\partial_jh_{ik}+\partial_kh_{ij}-\partial_ih_{jk}),\quad R^{[1]}=\partial_i\partial_jh_{ij}-\partial^2h_{ii}.
$$

√gRを部分積分したΓΓ形式の二次係数が式(W-01)の空間項となる。lapseの一次摂動はR^[1]と掛け合わさる。Fourier方向をzに取ればR^[1]=k²(h₁₁+h₂₂)、shiftはvを通して入り、欠落するconstraintはない。momentum方程式はp₁₃=p₂₃=p₃₃=0、Hamiltonianの線形条件はh₁₁+h₂₂=0に比例する。これらの四secondaryと四primaryは互いにPoisson可換である。

k≠0でgauge方向を商に取り、h₁₁=−h₂₂=A、h₁₂=Bとすると式(W-02)が従う。任意方向の非零格子波数では局所の直交基底で同じ代数となる。基底回転をtorusの格子自己同型と仮定する必要はなく、元のwavevectorと偏極tensorをそのまま戻す。

式(WD-02)。

$$
\begin{aligned}R+K^2-K_{ij}K^{ij}&=0,\qquad D_j(K^{ij}-g^{ij}K)=0,\\
K^{[1]}&=hI+\Sigma+\tfrac12\dot\gamma,\qquad\left\langle R^{[2]}\right\rangle=-\tfrac1{4q_0^2}\left\langle(\partial_a\gamma_{ij})^2\right\rangle.\end{aligned}
$$

zero-mean TTと空間一定のh、Σは空間平均で交差しない。二次補正metricに作用するR^[1]は全発散なので積分が零で、背景K=0では二次補正Kも二次拘束へ線形には入らない。従ってHamiltonian拘束の平均は6h²−trΣ²−〈dotγ²〉/4−〈(∂γ)²〉/(4q₀²)=0となる。

式(WD-03)。

$$
\mathcal J[\xi]=\int_\Sigma\pi^{ij}\mathcal L_\xi g_{ij},\quad \xi=\partial_a:\qquad \mathcal J^{[2]}[\partial_a]\ \propto\ \left\langle\dot\gamma_{ij}\partial_a\gamma_{ij}\right\rangle.
$$

背景π=0、g一定なので二次補正fieldの線形寄与が消え、平均momentumの残りは式(WD-03)となる。一定のmean metric・shearとのcross termは周期性で零である。これが式(W-03)の三つのmomentum条件である。h、Σは斉次線形方程式で一定であり、wave方程式と周期積分からenergy・momentumも一定となる。

この積分条件はコンパクト背景のlinearization-stability条件である。正のTT energyは純粋波動方向を排除するが、体積速度との混合方向を許す。続いて固定flat torusでの十分性と、明示的な二次補正を示す。

## 固定flat torusへの既知定理の適用

式(WD-04)。

$$
\partial_i\partial_j\nu-\delta_{ij}\Delta\nu=0,\qquad
\partial_iX_j+\partial_jX_i=0,\qquad
\det(\mathcal S_k^\dagger\mathcal S_k)=32(k^2)^5\quad(k\ne0).
$$

q₀=1の単位では随伴kernelはWD-04で与えられる。traceと周期性からνは定数、周期Killing方程式からXも定数であり、KIDは四本で尽きる。回転・boostの座標一次関数は固定格子に降下しない。最後の式は両symbolを合わせた四列行列についての非零mode対照である。

式(WD-05)。

$$
\begin{aligned}
\delta g&=C+2t(hI+\Sigma)+\gamma, &\delta\pi&=\Sigma-2hI+\tfrac12\dot\gamma,\\
\beta&=\mathrm{tr}(\delta\pi)=-6h,&\omega^{\rm TT}&=\Sigma+\tfrac12\dot\gamma,\\
\langle(\omega^{\rm TT})^2\rangle-\beta^2/6
+\tfrac14\langle(\partial\gamma)^2\rangle
&=\mathrm{tr}\Sigma^2-6h^2+\tfrac14\langle\dot\gamma^2+(\partial\gamma)^2\rangle.
\end{aligned}
$$

CとΣは一定、Σはtracefree、γは平均零TTである。WD-05はArms–Marsden–MoncriefのHamiltonian錐（同論文の式(3.5)）をW-03のenergy条件へ移す。三つのmomentum錐（同論文の式(3.1)）もWD-03と一致する。同論文の第3節は平坦な $`g_0`$ と $`\pi_0=0`$ の場合そのものを扱い、序論はKilling場4本の平坦時空 $`T^3\times\mathbb R`$ を含むと明記している。同論文第3節のTheorem 3.3（p.96）を、境界のないコンパクトCMC初期面の滑らかな周期データに適用する。例えば整数s≥5のHˢ×Hˢ⁻¹で曲線を解釈でき、最低正則性は主張しない。これにより四条件を満たす線形データは厳密拘束データ曲線の接方向となる。局所Cauchy発展とLC持ち上げからEC真空解族を得る。

これは同じ局所時間区間での初期値問題の結果であり、二時刻境界値問題、全時間正則性、有限Fourier空間だけで閉じる厳密族を意味しない。高次補正には新しいharmonicを許す。

## 明示二次補正と二波の閉包

W-04の単一定在波について、指数をεA₁+ε²A₂等と展開すると全10 Ricci成分の一次・二次係数が零となる。hの両符号で成立する。A₂,B₂,P₂を式に代入する局所検算であり、有限振幅で残差が厳密零という主張ではない。

式(WD-06)。

$$
\begin{aligned}
f&=a\cos x\cos t,\quad g=b\cos y\cos t,\quad h^2=(a^2+b^2)/24,\\
H_{\mu\nu}&=\mathrm{diag}(-6ht,\thinspace 2ht+g,\thinspace 2ht+f,\thinspace 2ht-f-g),\\
\mathcal K&=\lbrace (0,0),(\pm1,0),(0,\pm1),(\pm2,0),(0,\pm2),(\pm1,\pm1)\rbrace .
\end{aligned}
$$

式(WD-07)。

$$
\begin{aligned}
\ddot{\bar Q}_{\mu\nu,n}+|n|^2\bar Q_{\mu\nu,n}&=-2G^{[2]}_{\mu\nu,n}[H],\\
-\dot{\bar Q}_{0\nu}+\partial_i\bar Q_{i\nu}&=0,\\
\psi&=-\tfrac7{16}(a^2\cos2x+b^2\cos2y)-ab\cos x\cos y,\\
\Delta\psi&=-2G^{[2]}_{00}[H]{}(0),\\
\bar Q_{ij}(0)&=\delta_{ij}\psi,\qquad \bar Q_{00}(0)=\bar Q_{0i}(0)=0,\\
\dot{\bar Q}_{ij}(0)=\dot{\bar Q}_{00}(0)&=0,\qquad \dot{\bar Q}_{0i}(0)=\partial_i\psi.
\end{aligned}
$$

g=η+εH+ε²Qとし、barはηによるtrace reverseを表す。source保存とWD-07の初期条件により二次harmonic条件が保存される。13波数を全て保持すれば、強制oscillatorの積分でQを構成でき、局所Einstein成分・実場条件が二次まで成立する。公開補足のtwo_waves.jsonはsourceと補正係数を与える。斉次波数を除いた切断では同じ構成にならない。

[論文の入口](README_ja.md) / [参考文献](REFERENCES.md)

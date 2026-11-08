# 付録G. 拘束作用素・還元・再構成

## 共通の構成式

式(RG-01)。

$$
\begin{aligned}
\mathcal H[N]&=\int N\left[\frac{\pi^{ij}\pi_{ij}-\pi^2/2}{\sqrt g}-\sqrt g\thinspace R\right],\\
\mathcal H[X]&=\int\pi^{ij}\mathcal L_Xg_{ij},\qquad
Cz=\delta(\mathcal H_\perp,\mathcal H_i),\qquad
H=\delta^2\mathcal H[N],\\
S^{[2]}&=\int dt\thinspace [p^\dagger\dot h-\tfrac12z^\dagger Hz-\alpha^\dagger Cz].
\end{aligned}
$$

複素sectorは共役sectorと組み合わせて実作用を取る。CのHamiltonian行はδN/Nをmultiplierとするため背景Nを含む。RG-01を固定frameで変分し、非座標frameのLie微分と構造定数を保持する。例えばδΓ=G⁻¹(∇h+∇h−∇h)/2、δRic=∇δΓ−∇δΓから空間potentialのHessianを得る。これが公開行列を再生成する定義である。

式(RG-02)。

$$
\begin{aligned}
CJC^\dagger&=0,\quad \dot C+CJH=BC,\quad
CE=GE=0,\quad E^\dagger\Omega E=\Omega_r,\\
H_r&=E^\dagger HE-E^\dagger\Omega\dot E,\quad
\dot z=JHz+JC^\dagger\alpha,\\
\alpha&=-(GJC^\dagger)^{-1}(\dot GE+GJHE)y,\qquad z=Ey.
\end{aligned}
$$

rankが一定の独立拘束行を選んでGを固定する。RG-02のEの時間微分は背景が時間依存するため必須であり、特定のgaugeで補正が零でも省略を前提にしない。独立な4次元Einstein変分を再構成場へ代入し、空間成分と拘束に比例する時間成分を対照する。公開補足の係数はこの導出の検算・展開である。

## Nil base非零sector

式(NG-01)。

$$
\begin{aligned}
Gz=0:&\quad h_{12}=h_{13}=0,\quad h_{22}/A+h_{33}/B_0=0,
\quad2A\delta\pi^{11}+2\pi^{11}_0h_{11}=0,\\
P&=h_{33}/B_0,\quad Q=h_{23}/B_0,\quad B_0=16u^2/(9A),\\
\Pi_P&=-A\delta\pi^{22}+B_0\delta\pi^{33}
+(A\pi^{22}_0+B_0\pi^{33}_0)P,\\
\Pi_Q&=B_0p_{23}+2B_0\pi^{22}_0Q.
\end{aligned}
$$

式(NG-02)。

$$
\begin{aligned}
P_{uu}+u^{-1}P_u+9(k^2+2\zeta)P-18ik\zeta Q&=0,\\
Q_{uu}+(u^{-1}+2\bar P_u)Q_u+9k^2Q+18ikP&=0,\\
\bar P=\log\frac{4u}{3(1+u^4)},\quad \zeta=e^{2\bar P},\quad R=4u/3,&\qquad k\ne0.
\end{aligned}
$$

wavevectorをbaseの一軸へ合わせた表示である。整数波数を持つ全base方向には対応する偏極回転を戻す。選択した拘束minorは−128ik³u⁵(1+u⁴)/9、gauge minorは256ik³u/81で、区間全体で非零である。NG-01をRG-02へ代入すると本文N-01の作用を得て、NG-02へ一致する。空間gradientの非対角成分はcos/sinを混ぜる。k²−2ζ≥1/9という区間内の正定値対照は時間に依存しない保存energyを意味しない。

## Nil central非零block

NF-02の基底でRG-01を展開する。各Lで存在するmultiplierだけを含め、独立momentum列を拘束から解き、それに共役なmetric列をgaugeで固定する。存在する拘束行のrankは全て最大である。b=√μ>0とすると、L=−1,0,1,2の明示minorは次の通りである。L=−2には拘束行がない。

式(NG-03)。

$$
\begin{aligned}
M_{-1}&=2b(1+u^4),\qquad M_0=-256b^3u^5/9,\\
M_1&=-\frac{512b^3u^3[9b^2u^2(1+u^4)^2+12u^4+4]}{81(1+u^4)},\\
M_2&=-2048b^3u^5(1+u^4)/9.
\end{aligned}
$$

一般L≥2では六metric成分と四multiplierが存在する。有限閉包により随伴kernelは同じcentral非零sectorのKIDになるが、固定商の唯一のKIDはcentral零のZである。したがって四拘束は独立で二複素対となる。この構造論証と、五つのLでの明示行列計算を区別する。低端の減少はn<0成分の欠如でありKIDによるrank欠損ではない。

式(NG-04)。

$$
H_{r,-2}=-2\begin{pmatrix}
\dfrac{81\mu^2A^2-144\mu u^2+64}{72uA^2}&\dfrac{2(u^4-1)}{uA}\\
[4pt]
\dfrac{2(u^4-1)}{uA}&\dfrac{9A^2}{2u}
\end{pmatrix}.
$$

NG-04は一対だけの低端の明示例である。時間依存gauge上のHamiltonianの符号から一般のghostや不安定性は判定しない。L=2では還元、gauge multiplier、全10 Einstein成分を独立に照合した。他のLではこの4次元対照までを一律に行ったとはしない。

## Solのfiber零Gowdy sector

式(SG-01)。

$$
\begin{aligned}
ds^2&=-N^2dt^2+e^{2a}(dx+\beta dt)^2
+R[e^P(dy+Qdz)^2+e^{-P}dz^2],\\
R&=t+\rho,\quad P=2x+p,\quad Q=e^{-2x}r,\quad e^{2a_0}=F,\quad k=2\pi n/L_{\rm Sol},\\
C_H&=\delta p_R+(t-1/(4t))\delta p_a
+(2-1/(2t^2)-2k^2)\rho+2tikp,\\
C_x&=-2ik\delta a+ik(-2t-1/(2t))\rho+2\delta p_P-ik\delta p_a.
\end{aligned}
$$

n≠0でρ=δpₐ=0を選ぶとgauge determinantは−2ikである。δa=δpP/(ik)、δpR=−2tikpを消去しSO-01の作用を得る。pとrの方程式はそれぞれp̈+t⁻¹ṗ+k²p=0とr̈+t⁻¹ṙ+(k²+4)r=0となる。補助metric変数ℓ=4tṗ/(ik)、ℓ̇=4iktpも同時に再構成する。

## Sol一般非零orbitの微分右逆

以下は全orbit上の微分恒等式であり有限切断の逆行列ではない。h,pをU=diag(1,1,1,i,i,1)で変換し、拘束行をdiag(1,1,−i,−i)で実係数化する。primeはξ微分である。C=(Cq,Cp)のmomentum部分をCp p=rと書く。

式(SG-02)。

$$
\begin{aligned}
r_0&=A_0p_{11}+p_{22}+p_{33},\quad A_0=F(4t^2-1)/(2t),\\
r_1&=-2Fp'_{11}+2t(p_{22}-p_{33})+F(Kp_{12}+Mp_{13}),\\
r_2&=-t(p'_{12}+p_{12}+2Kp_{22}+Mp_{23}),\\
r_3&=-t(p'_{13}-p_{13}+2Mp_{33}+Kp_{23}).
\end{aligned}
$$

式(SG-03)。

$$
\begin{aligned}
f&=p_{11},\quad b=Kp_{12}-Mp_{13},\quad V=K^2+M^2,\quad d_0=2-FV/(2t),\\
a&=\frac{-b'-(K^2-M^2)r_0+A_0(K^2-M^2)f-Vr_1/(2t)-FVf'/t-Kr_2/t+Mr_3/t}{d_0},\\
p_{22}&=(r_0-A_0f)/2+r_1/(4t)+Ff'/(2t)-Fa/(4t),\\
p_{33}&=(r_0-A_0f)/2-r_1/(4t)-Ff'/(2t)+Fa/(4t),\\
p_{12}&=(a+b)/(2K),\quad p_{13}=(a-b)/(2M),\\
p_{23}&=(-r_2/t-p'_{12}-p_{12}-2Kp_{22})/M.
\end{aligned}
$$

SF-01とF/t≥e、V≥2|k₀l₀|からd₀≤2−4π²e/√5<0となり、区間1≤t≤2の全非零格子orbitで分母は消えない。r=0の写像をEp、f=b=0の写像をRpと定義する。TQ=(Q₁,0,0,KQ₂,−MQ₂,0)と置くと次の完全分解を得る。

式(SG-04)。

$$
\begin{aligned}
C_pE_p&=0,\quad C_pR_p=I,\quad T^\dagger E_p=I,\quad T^\dagger R_p=0,\\
E_pT^\dagger+R_pC_p&=I,\qquad Gz=R_p^\dagger h=0,\\
E&=\begin{pmatrix}T&0\\
-R_pC_qT&E_p\end{pmatrix},\quad
Q=E_p^\dagger h,\quad\Pi=(p_{11},Kp_{12}-Mp_{13}),\\
GJC^\dagger&=I,\quad GJG^\dagger=0,\quad E^\dagger\Omega E=\Omega_r.
\end{aligned}
$$

このgaugeはh₂₂=h₃₃=h₂₃=0、h₁₂/K+h₁₃/M=0と同値である。E†ΩĖはこの選択では厳密零だが、再構成にはĖが残る。RG-02からの再構成差をRとすると、拘束保存とgauge保存からCR=GR=0、還元方程式からE†ΩR=0となる。SG-04の完全分解でRはrange Eにあり、Ωrの非退化性によりR=0である。これは展開した全行列の総当たり相殺とは独立の証明である。

式(SG-05)。

$$
H_r=\sum_{a=0}^4 A_a(t,\xi)\partial_\xi^a,
\qquad H_r^\dagger=H_r,\qquad\mathrm{rank}\thinspace A_4=1,\qquad(J_rA_4)^2=0.
$$

公開補足はSG-05の全係数を展開し、最高階係数のrankと二乗零を確認する。これだけで双曲性・安定性を結論しない。

## S³低ラベルのrankと再構成

固定左添字あたりの拘束rankはj=0,1/2,1,3/2で3,8,11,16である。τ∈[1/4,1/2]における非零minorの分子の根数を有理数算術で数え、端点と分母も確認した。左重複度と実場条件を戻すと本文の完全block計数になる。j=0のℓ,m変分の二形式は[[0,−64],[64,0]]であり、時間変分はpure gaugeである。

式(TG-01)。

$$
\begin{aligned}
T_q&=\mathrm{diag}(1,1,1,i,1,i)\otimes\mathrm{diag}(1,\sqrt3,\sqrt3,1),\\
Gz=0:&\quad h_{11}=h_{22}=h_{33}=h_{13}=0\quad\text{for all right indices},\\
\det M&=\frac{2^{52}\tau^4(\tau-1)^4(\tau+1)^4
(\tau^4-6\tau^2+33)^2(3\tau^4+14\tau^2+3)^2}{81(1+\tau^2)^4}.
\end{aligned}
$$

TG-01はj=3/2の定数基底変換と使用したminorである。momentumは逆随伴で変換する。w=(1,√3,√3,1)として自由metricはh₁₂/(iw),h₂₃/(iw)、momentumは−iw p₁₂,−iw p₂₃の八対である。RG-02でEとmultiplierを構成し、独立4次元変分でδGij=0、δG00=−N CH/(2√g)、δG0i=−N Ci/(2√g)を確認する。全左添字を含めれば32実対となる。

## KIDと二次の必要条件

式(KG-01)。

$$
\int_\Sigma\left(|D\nu|^2+K_{ij}K^{ij}\nu^2\right)dV=0
\quad\Longrightarrow\quad\nu=0.
$$

選んだ各背景はCMCで $`K_{ij}K^{ij}>0`$ のため、KID lapseにKG-01を適用できる。残る空間Killing場は固定商に降下するものだけで、NilはZ一本、Solは零、S³ TaubはSU(2)×U(1)の四本である。Nilのcentral非零chargeはK-01となり、共役−mは同符号で加算する。全blockの和に条件を課し、block間の相殺を排除しない。m=0ではZ作用が零である。

式(KG-02)。

$$
L_a=-D_a^T,\qquad B_3=I_6\otimes D_3+S_3\otimes I_d,\qquad
(S_3h)_{ij}=-C^l{}_{3i}h_{lj}-C^l{}_{3j}h_{il}.
$$

S³の左三chargeはReΣ p̄ Lₐh、軸chargeはReΣ p̄ B₃hである。右scalar微分だけではcoframeのspinを落とすため、KG-02の両項が必要となる。これらはK-02から得る二次必要条件である。

式(KG-03)。

$$
\begin{aligned}
\left\langle2\delta p_P+\varepsilon\left[\delta p_P p_x+\delta p_r(r_x-2r)\right]\right\rangle&=0,\\
\langle\delta p_P^{(1)}\rangle=0,\qquad
2\langle\delta p_P^{(2)}\rangle+\left\langle\delta p_P^{(1)}p_x+\delta p_r(r_x-2r)\right\rangle&=0.
\end{aligned}
$$

Sol Gowdyの周期metric再構成では、 $`\Pi_PP_x+\Pi_QQ_x`$ のx平均 $`\langle\cdot\rangle`$ が零でなければならない（ $`\Pi_P,\Pi_Q`$ はSG-01の $`P,Q`$ の運動量）。 $`p,r`$ とその運動量 $`\delta p_P,\delta p_r`$ を $`\varepsilon`$ の一次とすると、この平均はKG-03の第1行の $`\varepsilon`$ 倍に等しい。第2行は $`\delta p_P=\delta p_P^{(1)}+\varepsilon\delta p_P^{(2)}+O(\varepsilon^2)`$ と展開したもので、斉次な二次補正を一次の場で定める。一般orbitのKID chargeではない。NilとS³ Taubでは、ここに示したcharge零だけから非線形持ち上げの十分性を結論しない。SolのKID消失に基づく十分性は、次の関数空間と定理適用により確立する。

## Solの局所可積分性

### 背景と全大域KID

SF-01の商は $`A\in SL(2,\mathbb Z)`$ によるtorus mapping torusであり、連結・コンパクト・境界なしの向き付き3次元多様体である。大域coframe $`(dx,e^xdy,e^{-x}dz)`$ 上のBF-01は全 $`t>0`$ で滑らかで、N>0、G正定値である。各 $`t_0\in[1,2]`$ を含む開時間区間の積時空では、空間のコンパクト性と因果速度の局所有界性により各sliceがCauchy面となる。付録Fの真空性と、次の外部曲率を用いる。

式(SL-02)。

$$
\begin{aligned}
K^i{}_j&=\frac1{\sqrt F}\mathrm{diag}
\left(t-\frac1{4t},\frac1{2t},\frac1{2t}\right),\\
\tau:=\mathrm{tr}\thinspace K&=\frac{t+3/(4t)}{\sqrt F},\qquad
|K|^2=\frac{(t-1/(4t))^2+1/(2t^2)}F>0.
\end{aligned}
$$

KIDを(f,X)とする。Kの正符号規約でmetric逆行列の変分も含めてtrace発展を取ると、真空拘束上で $`0=\Delta f-|K|^2f+X(\tau)`$ となる。ここで $`\Delta=\mathrm{div}\nabla`$ 。背景はCMCなのでX(τ)=0であり、境界なしの部分積分がKG-01を与え、f=0となる。残る条件は $`\mathcal L_Xg=\mathcal L_XK=0`$ である。

空間Ricci endomorphismの固有値は(−2/F,0,0)で、負の固有方向は単純である。Killing流はこの方向の単位vectorを保つため $`[X,\partial_x]=0`$ 。残るKilling方程式を解くと、普遍被覆の全局所Killing場は $`a(\partial_x-y\partial_y+z\partial_z)+b\partial_y+c\partial_z`$ と尽くされる。二つのfiber格子並進への降下条件からa=0、base deck変換から $`(b,c)=(\lambda^{-1}b,\lambda c)`$ となりb=c=0である。従って大域KIDは零となる。離散等長変換の不存在は必要としない。

### 関数空間と全射性

整数s≥6を固定し、Densを実density bundleとして次の空間を取る。Sobolev normと形式随伴の内積は滑らかな背景g₀で固定する。

式(SL-03)。

$$
\begin{aligned}
\mathcal X_s&=H^s(S^2T^{\ast}\Sigma)\times
H^{s-1}(S^2T\Sigma\otimes\mathrm{Dens}),\\
\mathcal Y_s&=H^{s-2}(\mathrm{Dens})\times
H^{s-2}(T^{\ast}\Sigma\otimes\mathrm{Dens}),\\
\Phi(g,\pi)&=(\mathcal H_\perp,\mathcal H_i),\qquad
P=D\Phi_{(g_0,\pi_0)}=(P_g,P_\pi).
\end{aligned}
$$

拘束写像はRG-01の全Hamiltonian・momentum拘束であり、正定値metricの開集合で $`\Phi:\mathcal X_s\to\mathcal Y_s`$ は滑らかである。densityの同定には固定背景の体積densityを用いる。 $`\pi^{ij}=\sqrt g(K^{ij}-g^{ij}\tau)`$ とその逆はこの正則性で滑らかであり、pはδπを表す。[FischerMarsdenMoncrief1980](REFERENCES.md#fischermarsdenmoncrief1980) p.155の空間設定を採用し、最低正則性は主張しない。

Hamiltonian行はgに二階・πに零階、momentum行はgに高々一階・πに一階である。momentum density変分のbundle $`S^2T\Sigma\otimes\mathrm{Dens}`$ に作用する正の二階作用素Bを用いて、次の空間楕円作用素を作る。

式(SL-04)。

$$
\begin{aligned}
B&=I+\nabla^{\ast}\nabla,\qquad
\mathscr L=P_gP_g^{\ast}+P_\pi BP_\pi^{\ast},\\
\sigma_4(\mathscr L)(\xi)&=2|\xi|^4
\mathrm{diag}\left(1,I_3+|\xi|^{-2}\xi\otimes\xi\right),\\
\langle w,\mathscr Lw\rangle&=\lVert P_g^{\ast}w\rVert^2+\lVert P_\pi^{\ast}w\rVert^2
+\lVert\nabla P_\pi^{\ast}w\rVert^2.
\end{aligned}
$$

normは固定背景のL² norm、 $`\xi\ne0`$ は空間covectorである。背景正規直交frameでHamiltonianの主symbol随伴は符号を除き $`f(\xi\otimes\xi-|\xi|^2I)`$ 、momentumのπ主symbolは $`p\mapsto-2ip(\xi,\cdot)`$ 。非対角blockは高々三階なので、SL-04の四階主symbolは正定値となる。滑らかなwで $`\ker\mathscr L=\ker P^*`$ 。distributionとして $`P^{\ast}w=0`$ なら $`\mathscr Lw=0`$ であり、楕円正則性によりwは滑らかになるため、上のKID消失が弱いkernelにも及ぶ。

閉多様体上の自己随伴楕円作用素のFredholm理論により $`\mathscr L`$ は可逆で、次が有界右逆と補空間を与える。

式(SL-05)。

$$
\begin{aligned}
\mathscr R&=\begin{pmatrix}P_g^{\ast}\\
BP_\pi^{\ast}\end{pmatrix}
\mathscr L^{-1}:\mathcal Y_s\longrightarrow\mathcal X_s,
\qquad P\mathscr R=I,\\
\mathcal X_s&=\ker P\oplus\mathrm{range}\thinspace\mathscr R,\qquad
\mathrm{proj}_{\ker P}=I-\mathscr RP.
\end{aligned}
$$

$`\mathscr L^{-1}`$ は $`H^{s-2}`$ を $`H^{s+2}`$ へ写し、右逆の二成分はそれぞれ高々二階・三階なので、値は正しく $`\mathcal X_s`$ へ入る。これにより全射性・像の閉性・kernelの分裂を得る。これは同文献Theorem 2.1の楕円的分解を本背景で具体化したものである。空間楕円補助作用素 $`\mathscr L`$ は、4階時間発展系SG-05のHamiltonianとは異なる。

### 厳密族・局所発展・LC持ち上げ

x₀=(g₀,π₀)とし、 $`u\in\ker P`$ 、 $`v\in\mathcal Y_s`$ に対して $`\Phi(x_0+u+\mathscr Rv)=0`$ を考える。v微分は恒等写像なので陰関数定理からv=ψ(u)、ψ(0)=Dψ(0)=0を得る。任意のz∈ker Pについて $`x_\varepsilon=x_0+\varepsilon z+\mathscr R\psi(\varepsilon z)`$ がSL-01を満たす。初期データ族はεについて滑らかである。CMCは背景でKIDを除く条件であり、近傍の全データに固定平均曲率を課していない。滑らかさを述べる対象は拘束零集合であって、離散対称性も含む全微分同相商ではない。

Harmonic/wave-map gaugeによる局所Cauchy理論は、十分小さいεの族に共通する短時間の真空発展を与える[HughesKatoMarsden1977](REFERENCES.md#hugheskatomarsden1977)。拘束とgaugeの伝播を含めて全真空方程式を満たし、一意性は微分同相を除いて述べる。発展は $`C^0H^s\cap C^1H^{s-1}`$ に属し、接方向は一階低い $`C^0H^{s-1}\cap C^1H^{s-2}`$ 位相で解釈できる。最高位Sobolev normで発展写像が滑らかであることまでは要求しない。近傍metricの平方根からcoframeを選びLC接続を付ければ、A-02により同じEC＋定数完全NYの局所解族となる。接続はcoframeの一階微分を含む。この初期値問題は二時刻の固定境界データを課す問題ではない。

### 個別modeの適用範囲

Gowdyの非零base modeでは、SG-01の $`\delta a=\delta p_P/(ik)`$ 、 $`\delta p_R=-2tikp`$ と全metric・momentumへの再構成が滑らかな商上の拘束データを与える。実場には反対波数の共役を戻す。有限和に加え、平均零の周期場p,r∈Hˢ、運動量δpP,δpr∈Hˢ⁻¹にも適用できる。逆微分がδaの正則性を一階上げ、δpRはpの一階微分となるためである。k=0にはこの割算を適用しない。KG-03の二次平均補正を許すこととSL-01は整合する。

一般orbitでは、SF-02の等長写像と大域frameを用いる。

式(SL-06)。

$$
\begin{aligned}
X_iU_qf&=U_qD_if,\qquad
U_q\mathcal D\subset\bigcap_{m\ge0}H^m(\Sigma)\subset C^\infty(\Sigma),\\
E(t_0)\mathcal D^4&\subset\mathcal D^{12},\qquad CE=0.
\end{aligned}
$$

各有限word $`D_{i_1}\cdots D_{i_m}`$ は指数関数の重み付き微分の有限和なので、 $`\mathcal D`$ の定義によりその結果はL²に入る。閉多様体上のframe微分normとSobolev normの同値性から第1行を得る。重み付きnormでcompact-support cutoff近似を取れば、弱微分の交換も正当化される。SG-03の分母は零から離れ、再構成係数とその各微分は高々指数増大なのでEは $`\mathcal D`$ を保つ。CE=0と元のtensor成分への位相復元・共役orbitの加算により、実初期データはSL-01の対象となる。

SH-03の有限Hermite和は多項式とGaussianからなり $`\mathcal D`$ に属するため、その実数直線全体の連続初期profileも含む。有限窓の離散軌道を厳密解と認定する結果ではない。複数orbitの和には、共役を処理した再構成成分について

式(SL-07)。

$$
\sum_{[q]}\left(\lVert h_{[q]}\rVert_{H^s(\Sigma)}^2
+\lVert p_{[q]}\rVert_{H^{s-1}(\Sigma)}^2\right)<\infty.
$$

を要する。異なるorbitの直交性とframe微分によるorbit保存から全体が $`\mathcal X_s`$ で収束し、Pの連続性により全拘束を満たす。各orbitが個別に滑らかであることや、還元係数の重みなし二乗和だけからSL-07は従わない。任意のL²還元データ、SG-05の閉作用素実現・well-posedness・数値収束、全時間安定性、単一sector内の非線形閉包は本結果に含めない。

[論文の入口](README_ja.md) / [参考文献](REFERENCES.md)

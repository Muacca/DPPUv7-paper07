# 3. 背景と六成分の斉次作用

式(H-01)。

$$
\begin{aligned}c&=\frac{\mathcal V_{\Sigma}}{2\kappa^2},\quad d=\det g,\quad v=\dot g-\mathcal L_{\beta}g,\quad K=\frac1{2N}g^{-1}v,\\
P(n,g)&=\frac12[\mathrm{tr}(ng)]^2-\mathrm{tr}(ngng),\qquad {}^{(3)}R=\frac{P(n,g)}{d},\\
L_0&=c\sqrt d\left[\frac{\mathrm{tr}((g^{-1}v)^2)-[\mathrm{tr}(g^{-1}v)]^2}{4N}+N\frac{P(n,g)}d\right].\end{aligned}
$$

独立接続方程式を解き、境界を処理した後の作用である。gの六成分、lapse、三shiftを変分してからgaugeを選ぶ。shiftは固定商へ降りる不変ベクトル場を用いる。class Aでは不変体積の発散項が消え、この全斉次変分はhomogeneous Einstein方程式を保持する[Torre2004](REFERENCES.md#torre2004)。対角計量だけに先に制限した計数は本節の計数ではない。

曲率式は、g=EᵀE、物理的構造行列m=EnEᵀ/det E、およびR=½(tr m)²−tr(m²)から得る。kinetic項はKᵢⱼKⁱʲ−K²である。一般の正定値gでLegendre写像は可逆で、六速度のDeWitt形式は五正・一負の符号を持つ。この符号だけから物理ghostを判定しない。

式(H-02)。

$$
\begin{aligned}\pi^{ij}&=c\sqrt d\thinspace (K^{ij}-(\mathrm{tr}K)g^{ij}),\qquad T(g,\pi)=\mathrm{tr}(\pi g\pi g)-\tfrac12[\mathrm{tr}(\pi g)]^2,\\
\mathcal H&=\frac{T-c^2P}{c\sqrt d},\qquad F=c\sqrt d\thinspace \mathcal H=T-c^2P.\end{aligned}
$$

πは空間積分済みの対称momentumで、πⁱʲδgᵢⱼは全添字和である。非対角成分を独立座標にすると、その共役運動量は2πⁱʲとなる。Fへの置換は正の非零因子c√dによる拘束の再表示であり、零集合上でHamiltonian方向を同じ非零因子だけ変える。

## 共通scale背景の判定

Paper06の共通scale背景の判定を、本稿の拘束解析の出発点として式(B-01)にまとめる[paper6](REFERENCES.md#paper6)。

式(B-01)。

$$
\begin{aligned}g&=q^2 I,\quad 6\chi=\tfrac12(\mathrm{tr}n)^2-\mathrm{tr}(n^2),\quad X=2n^2-(\mathrm{tr}n)n+4\chi I,\\
G_{00}&=3(H^2+\chi/q^2),\qquad G_{0i}=0,\\
\tfrac13\sum_iG_{ii}&=-(2DH+3H^2+\chi/q^2),\qquad (G_{ij})_{\rm TF}=X_{ij}/q^2.\end{aligned}
$$

全接続shell上では真空Einstein方程式となる。共通scale縮約の実軌道は(Dq)²+χ=0を満たすが、NilとSolではXが残る。標準NilはX=diag(−1,−1,2)/3、Solはdiag(−4,2,2)/3である。静的torusは完全解であり、round sphereにはこのモデルの実共通scale真空軌道がない。

したがって静的sphereの二次作用はoff-shellな代数的対照であり、真空粒子のスペクトルではない。Nil・Solのshape固定軌道にも同じ注意が必要である。一般の異方的真空解の存在を否定しているのではない。第5節の正則領域は一般gとmomentumを用い、付録Cに拘束を満たす代表点を示す。

摂動はe=ē+εf、ω=ω̄+εaから同じ次数で展開する。一次変分を捨ててon-shellとみなさず、Noether恒等式の背景残差を保持する。Nil・Solの等しいcoframe scaleは空間の等方対称性ではなく、scalar/vector/tensor blockの独立閉包を保証しない。

## 摂動に用いる時間依存真空背景

曲がった幾何では共通scale参照計量を用いず、次の時間依存真空背景を採用する。Gは付録Aの固定不変coframeでの空間metricである。全4次元Ricci、ADM拘束、固定商への降下を確認する構成と規格化を付録Fに示す。

式(BF-01)。

$$
\begin{aligned}
\mathrm{Nil}:&\quad G=\mathrm{diag}(A,A,16u^2/(9A)),\quad N=3\sqrt A,\quad A=1+u^4,\\
\mathrm{Sol}:&\quad G=\mathrm{diag}(F,t,t),\quad N=\sqrt F,\quad F=e^{t^2}/\sqrt t,\\
S^3:&\quad G=\mathrm{diag}(4(1+\tau^2),4(1+\tau^2),16(1-\tau^2)/(1+\tau^2)),\\
&\quad N=\sqrt{(1+\tau^2)/(1-\tau^2)}.
\end{aligned}
$$

NilとSolはそれぞれ1≤u≤2、1≤t≤2、S³ Taubは1/4≤τ≤1/2を用いる。全区間でN>0、G正定値である。S³ではℓ=1、Taub parameter m=0を選び、Cauchy horizonと時間対称点を避ける。これらは全接続方程式でLC接続へ持ち上げた同じEC真空背景である。

[論文の入口](README_ja.md) / [参考文献](REFERENCES.md)

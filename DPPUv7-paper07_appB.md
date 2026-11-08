# 付録B. 全contortion・Euler式・端点

式(CB-01)。

$$
\begin{aligned}K_{ab}&=K_{ab,c}e^c,\qquad K_{ab,c}=-K_{ba,c},\qquad s_a=\eta_{aa},\\
R_{\rm EC}&=R_{\rm LC}+\nabla_aJ^a+R_K,\qquad J^a=2s_a\sum_bs_bK_{ab,b},\\
R_K&=\sum_{a,b,e}s_as_bs_e(K_{ae,a}K_{eb,b}-K_{ae,b}K_{eb,a}).\end{aligned}
$$

a<b、c=0,…,3を独立な24座標として固定する。二次形式R_KのHessianはdet=2³²、rank=24で、KからTへの24×24線形写像も可逆である。このdetは座標規格化に依存するが、非退化性は依存しない。第2節のCartan方程式の内部積による証明が一般の非退化coframeでの一意性を与え、この成分計算は独立な代数対照となる。

shiftを含む一般のcoframeでも式(CB-01)はtensor恒等式である。spatial sliceに境界はなく、発散の積分は時間端点の法線currentのみとなる。LC GHYとこのcurrent、定数NYを対応する境界代表で処理すると、S₀=S_ADM+∫√−g R_K/(2κ²)を得る。K方程式は点ごとにK=0を与える。

式(CB-02)。

$$
\begin{aligned}\omega&=K+W(e,\partial e),\quad W=\omega_{\rm LC},\\
E_K(S')&=E_\omega(S),\\
E_e(S')&=E_e(S)+E_\omega(S)\frac{\partial W}{\partial e}-\partial_\mu\left(E_\omega(S)\frac{\partial W}{\partial(\partial_\mu e)}\right).\end{aligned}
$$

添字の縮約を省略して書いたEuler chain ruleである。完全な接続方程式を課せば追加項が消えるため、coframe方程式もLC Einstein方程式へ対応する。接続shellを使う前にEuler式の一部だけを捨てることは正当化されない。

変分の端点1形式には、式(CB-02)の部分積分によるE_ω比例項が加わる。それは接続shellへの引戻しで消える。さらに境界代表の差はδB_endであり、δ²B_end=0から還元二形式が一致する。これは微分依存の変数変換をpoint正準変換とみなす主張ではなく、解空間への引戻しと対応する端点データについての主張である。

式(CB-03)。

$$
\begin{aligned}L_0&=c\left[-\frac{6q\dot q^2}{N}+6\chi Nq-6Nq\eta^2+\frac23Nq^3V^2\right],\\
B_{\rm end}&=-2cq^3V+6\theta\mathcal V_{\Sigma}q^2\eta,\\
(p_\eta,p_V)&=(0,0),\qquad(\eta,V)=(0,0).\end{aligned}
$$

共通scale族では最後の四条件が四つの第二級拘束となる。端点を保持する代表ではprimaryはP_η−∂_ηB=0、P_V−∂_VB=0へ、q momentumもP_q=p_q+∂_qBへ移る。拘束とHamiltonianを同じ置換で移せばPoisson代数・Dirac二形式は変わらない。θが定数でなくなる場合やspin source、退化coframeにはこの結論を拡張しない。

## 非斉次blockでの直接対照

式(CB-04)。

$$
B_{\rm EH}=\tfrac12\epsilon_{abcd}e^a\wedge e^b\wedge K^{cd},
\qquad B_{\rm NY}=e^a\wedge T_a,\qquad
dB_{\rm NY}=T^a\wedge T_a-e^a\wedge e^b\wedge R_{ab}.
$$

式(CB-04)では結合係数を除いたcurrentを表示した。EH currentはcontortion座標、NY currentは生の接続摂動で展開する場合があるため、両者を足す前に同じ接続へ代入する。二次摂動ではcos/sin間の積と、背景との積を通じる必要な平均二次Cartan項を保持する。

flat非零modeではLorentz六方向を含む40配置変数の連鎖、Lorentz固定34変数、計量10変数を直接比較した。Nil baseとSol Gowdyでは全24成分Cartan、補助二次形式、端点currentを直接検算した。S³のj=3/2では全右添字のCartanと二次端点恒等式を確認した。Nil centralとSol orbitではtensor/orbit閉包を保つCartan解を確認したが、追加の全Dirac連鎖・全二次端点行列を計算した範囲とは区別する。

[論文の入口](README_ja.md) / [参考文献](REFERENCES.md)

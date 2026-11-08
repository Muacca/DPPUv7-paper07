# 付録C. 拘束rankとroute対照の詳細

一定の対称nはSO(3)で対角化できる。式(H-03)の行列Aは構造代数のadjoint作用であり、Jacobi恒等式からmomentumのPoisson代数が従う。TとPはこの作用で不変なので{J,F}=0となる。六座標とその共役を用いれば、これらは多項式恒等式として直接検算できる。

式(CC-01)。

$$
\begin{aligned}g&=\mathrm{diag}(1,2,3),\quad c=1,\quad M=\pi g,\\
M_{T^3}&=\mathrm{diag}(1,1,0),\quad M_{\rm Nil}=\sqrt3 I,\\
M_{\rm Sol}&=\frac{5\sqrt3}{3}I,\quad M_{S^3}=\mathrm{diag}(4\sqrt2,0,0).\end{aligned}
$$

全てF=J=0を満たす。T³のminorはdF/dp₃₃=−6である。他の代表では行(F,J₁,J₂,J₃)の必要な先頭部分、列(p₁₁,p₂₃,p₁₃,p₁₂)の必要な先頭部分を取ると、Nilでは9√3、Solでは50√3、S³では−64√2となる。Solの値は対角化した $`n=\mathrm{diag}(0,1,-1)`$ で評価する。従って拘束Jacobianのrankは1/3/4/4である。

Aの線形包の次元は0/2/3/3で、Jのrankの上限を与える。非零minorの近傍ではこの最大rankとdFの独立性が保持され、式(H-04)の正則層となる。Σの等方化や対称性増大などでminorが消える点へ一律に計数を延長しない。

## 二次作用のclosure

式(CC-02)。

$$
M_{A b}=\lbrace \Phi_A,\phi_b\rbrace ,\quad u^AM_{Ab}=0\quad\Longrightarrow\quad u^A\lbrace \Phi_A,H_c\rbrace \simeq0.
$$

φは元のprimary、Φはその段階までの全独立拘束である。Mの左nullspaceから新しい独立線形拘束だけを追加し、追加がなくなるまで反復する。定数係数二次作用ではこの線形代数手順がexactである。非線形拘束を線形係数rankへ落とす用途には使わない。最終Poisson行列rankをs、独立拘束総数をmとすれば第一級はm−s、第二級はsとなる。

| 対照 | 接続：第一級/第二級/pair | 消去後：第一級/第二級/pair |
| --- | --- | --- |
| 完全作用torus二次 | 1 / 4 / 1 | 1 / 0 / 1 |
| 完全作用sphere二次 | 0 / 8 / 0 | 0 / 4 / 0 |

sphereの消去後はH=−p_u²/24−6ℓu、p_ℓ=0からu=0、p_u=0、ℓ=0を得る。torusではlapse primaryが第一級となり、二次作用にはscalar pairが残る。元の非線形拘束による条件は第5節で扱う。

## 完全作用からのscalar対照

c=1、n=2kIとし、独立接続をω⁰ᵢ=b eⁱ、ωʲₖ=a eⁱで表示する。k=0はtorus、k=1はround sphereの参照である。

式(D-02)。

$$
L_{\rm P}=6Nq^3b^2-12q^2\dot q\thinspace b+12kNq^2a-6Nq^3a^2,\quad a=k/q,\quad b=\dot q/(Nq).
$$

式(D-03)。

$$
\begin{aligned}L_{\rm met}&=6k^2Nq-6q\dot q^2/N,\\
q&=1+\varepsilon u,\quad N=1+\varepsilon\ell,\quad a=k+\varepsilon A,\quad b=\varepsilon B,\\
L^{[1]}&=6k^2(\ell+u),\\
L_{\rm P}^{[2]}&=-6A^2-12kAu+6B^2-12B\dot u+6k^2\ell u-6k^2u^2,\\
L_{\rm met}^{[2]}&=-6\dot u^2+6k^2\ell u.\end{aligned}
$$

二次接続方程式A=−ku、B=dot uを代入すると同じmetric作用になる。sphereでは一次変分が非零で、両経路のpair数0はoff-shellな代数的対照である。torusの二次scalar作用には1対が残るが、元のlapse方程式はdot u²=0を要求する。これは第5節の混合錐のpure-scale切断に対応する。

## 対称Taub点と非零modeの計数

BF-01のS³ Taub点ではmomentum rank 2、全rank 3であり、斉次線形商は3対である。正則代表CC-01のrank 3、4とは別の層である。j=0とj=1のKIDによる例外、および区間全体の非零minorは付録Gに示す。flat非零modeの40/34/10配置変数の比較は第4節の一つの実Fourier成分を単位とし、斉次正則層のH-04とは別の計数である。

[論文の入口](README_ja.md) / [参考文献](REFERENCES.md)

# 規約と以前のlocal blockとの対応

現在の作用・固定商・物理半径は本文第2節・付録A–Bを用いる。旧torsion入力配列は幾何学的torsionの負号 $`\tau=-T`$ 、NYの旧記号は $`\theta=\theta_{\mathrm{NY}}/\kappa^2`$ に対応する。構造係数の大きさがa/rの球面ではq=2r/|a|であり、a=1、2、4を同時に換算する。半径・体積・波数の単位は一組として扱う。

付録Iの式(D-04)のtrace blockは、coframe摂動の対角成分をs/3、lapse摂動をℓとし、接続摂動の01,02,03対角成分をb、12_3と23_1をa、13_2を−a、その他を零としたものに対応する。境界項を合わせて−8s dot bを8 dot s bへ移すと、a=−s/6−ℓ/2、b=−dot s/3を消去でき、さらにℓ=s/3からL=−4 dot s²/3となる。

このstrict二次blockは、接続表示で0第一級・6第二級・1 pair、消去後で0第一級・2第二級・1 pairとなる。元の26配置変数blockは0第一級・50第二級・1 pairである。primary保存では全てのprimary-Hamiltonian bracketを零と置かず、primary Poisson行列のmultiplier項を含める。旧local actionに含まれたlapse摂動×接続二乗は三次項として分ける。

このblockと付録Cの完全作用由来のscalar二次作用は区別する。`routes.json` の `legacy_trace` はこの比較、`0`・`1` は完全作用からのtorus・sphere対照に対応する。付録D・F–Hも同じ作用規約に従う。版間の科学的変更は[CHANGELOG](../CHANGELOG.md)、現在の再現手順は[script guide](README.md)を参照する。

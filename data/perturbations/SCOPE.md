# Adopted results and unresolved scope / 採用結果と未確立の範囲

This table accompanies v2. 「解析」「既知定理の適用」「有限時間数値例」を同じ証拠として扱わない。

| Result | Adopted statement / 採用範囲 | Boundary / 未確立・適用外 |
| --- | --- | --- |
| Flat torus, four Taub conditions | Necessary and sufficient for local integrability under Appendix D's compact CMC and regularity hypotheses; an application of AMM Theorem 3.3 | Not proved by finite symbolic checks; no all-time or fixed-finite-harmonic exact-family claim |
| Single and nonparallel waves | Explicit local vacuum corrections through order epsilon², with all 13 generated wavevectors for two waves | Not exact finite-amplitude solutions |
| Independent connection | Flat nonzero real quadrature: 40/34/10 configurations, 20/8/8 first-class and 36/48/0 second-class constraints, two pairs | Direct full flat-chain result; curved Cartan and endpoint checks have the narrower scopes of Appendix B |
| Regular homogeneous metrics | T³/Nil/Sol/S³ = 5/3/2/2 pairs | Static torus and symmetric Taub points require separate rank analysis |
| Nil base | Full reduced action and reconstruction on the selected LRS vacuum; complete real opposite-wavevector block has four pairs | No non-LRS or nonlinear lifting sufficiency result |
| Nil central | Finite tensor-grade closure; one complex pair at L=−2,−1,0,1 and two at L≥2 per positive m and lattice copy | Five explicit grades; independent 4D reconstruction at L=2; lower-end missing slots are not KID defects |
| Sol Gowdy | Fiber-zero, base-nonzero two-field action and periodic reconstruction | The k²+4 coefficient is not a conserved physical frequency; periodic reconstruction includes a homogeneous second-order correction |
| Sol general orbit | Exact infinite-orbit closure, differential right inverse and two complex canonical field pairs on the stated common core | No finite Hermite closure, closed-operator extension or whole-line wellposedness theorem for the reduced system SG-05 |
| Taub Wigner blocks | Real ranks 3/16/33/64 and pairs 3/8/21/32 for j=0,1/2,1,3/2; complete j=3/2 reconstruction | j=0,1 have KID defects; no all-j or general triaxial classification |
| Nil and Taub KIDs | Nil 1 and Taub 4; explicit quadratic necessary charges | Nonlinear lifting sufficiency is not established |
| Sol local integrability (`PROVED_SCOPED`) | No global KIDs on the fixed BF-01 quotient; every real full linearized constraint datum in H^s x H^(s-1), integer s >= 6, is tangent to an exact vacuum constraint family with local Cauchy development and LC lift (SL-01–05) | Small family parameter and local time; no fixed CMC value, gauge, Gowdy ansatz or single-orbit nonlinear closure is imposed; no uniform guarantee on [1,2] or nonlinear stability |
| Sol mode applicability | Gowdy nonzero-base reconstruction; real general-orbit data E(t0)y with y in D^4, original tensor phases and conjugate orbit restored; the continuous finite Hermite profile SH-03; infinite sums satisfying SL-07 | Arbitrary reduced L2 data and finite-window numerical trajectories are not certified; the theorem is distinct from closed realizations or whole-line well-posedness of SG-05 |
| Finite blocks | Bessel, second-order/Hamilton and reduced/full-ADM comparisons for the stated initial data | Finite-time linear examples, not exact curved nonlinear families or stability theorems |
| Sol Lorenz example | One orbit (1,0), one initial profile, t in [1,2], distance at most 4 from the center; separate space/window/time comparisons and actual dense-polynomial residuals | Normalized strong residual (SH-05) of the reduced fourth-order system is 1.54e−3 in the dense test and does not decrease under spatial refinement; no convergence in that form, whole-line guarantee or physical spectrum |

Solでは、全拘束データの解析的な局所可積分性（SL-01–07、`PROVED_SCOPED`）と、有限時間・内部領域の独立Lorenz数値対照を別の根拠として採用する。4階還元系での収束・閉作用素実現、任意のL²還元データのwell-posedness、Nil・Taubの非線形持ち上げ十分性、全orbit・全j・全時間安定性、固定sector内の厳密族、非LRS Nil、一般三軸Bianchi IX、特異拘束層の全分類は未確立である。固定frame成分のL²残差と、高階還元式の正規化強残差は別のnormであり、数値を直接比較して同じ判定にしない。

The mathematical derivations are in Sections 4–7 and Appendices D, F–H, not only in code. Geometric Laplacians in Appendix E remain separately labeled reference operators.

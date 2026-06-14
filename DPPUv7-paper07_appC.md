# DPPUv7-paper07 — Appendix C: Route Comparison and Route-Dependent Discrepancy

[← Appendix B](DPPUv7-paper07_appB.md) | [Appendix D →](DPPUv7-paper07_appD.md)

---

## Appendix C. RBD route comparison and route-dependent discrepancy

本 appendix では [Section 3](DPPUv7-paper07_sec03.md) の対照実験の詳細を示す。ここでの "Palatini reduction non-commutativity" は route-level の記述用語であり、一般定理の主張ではない。

### C.1 The RBD route

RBD（Routh-before-Dirac）route は、空間接続の補助成分（spatial-omega）をその場の方程式で解いて作用に代入してから、残った変数で制約解析を行う手続きである。この代入は auxiliary variable の条件（微分を含まない代数的可解性）の検証を経ておらず、一般論 [13] が警告する非可換性のリスクを負う。

### C.2 Direct route recomputation

本稿の negative control は toy model ではなく、同梱の DBR/RBD 制約力学エンジン（`dppu.dirac`）を用い、同一の $S^3$ Branch-A local first-order Lagrangian から route variants を再計算するものである。再現スクリプト `script/scripts/paper07/route_comparison_s3_branch_a_witness.py` は中間 CSV を証拠として読まず、以下を直接実行する。

- **DBR_KEEP_ALL**: 独立接続を保持したまま Dirac 閉包を行う。
- **RBD_SPATIAL_OMEGA_ROUTH**: spatial-omega を先に Routh 消去し、その後に scalar projection と Dirac 閉包を行う。
- **RBD_FULL_K_ROUTH**: full- $K$ を先に Routh 消去し、その後に scalar projection と Dirac 閉包を行う。

結果は次の通りである。

| route | $N_{\rm phys}$ at $l=0/1/\ge2$ | $l=0$ surviving pair | $H_{PP}$ at $l=0$ | status |
|---|---:|---|---:|---|
| DBR_KEEP_ALL | $0/0/0$ | none | n/a | computed |
| RBD_SPATIAL_OMEGA_ROUTH | $1/1/2$ | $(s_{l0},P_{s,l0})$ | $-3/8$ | computed |
| RBD_FULL_K_ROUTH | $1/1/2$ | $(s_{l0},P_{s,l0})$ | $-3/8$ | computed |

したがって route 非可換性の主証拠は $N_{\rm phys}$ の不一致であり、 $l=0$ の $H_{PP}=-3/8$ は RBD 側の余分 scalar branch が benign な再表現でないことを示す診断である。ただし、RBD は採用しない route なので、この負符号を physical ghost / physical instability の主張には用いない。

### C.3 Comparison setup and result

対照実験は $S^3$ 背景・Branch-A scalar sector で実施した。§C.2 の直接再計算は、この sector の操作順依存性を DBR/RBD エンジン上で再現するものである。

- **DBR 側**: Dirac 閉包後、物理的 scalar 自由度は全 harmonic レベルで 0（kernel は $S^3$ では不在; $d_{\rm zero}=0$ ）。
- **RBD 側**: spatial-omega Routh と full- $K$ Routh はともに $N_{\rm phys}=1/1/2$ を返し、DBR と非可換である。
- **同一 Lagrangian 内部実証**: $l=0$ では同一の Branch-A local first-order Lagrangian から出発して、DBR: 0 / RBD: 1 の不一致と $H_{PP}=-3/8$ が厳密に（数値誤差なしの symbolic 計算で）再現される。これにより不一致が実装差や正規化差でなく route 自体に起因することが確定する。
- **スコープ限定**: $l\ge1$ の RBD extra count は、Routh 先行が Lorentz gauge 構造を壊す route-level artifact として記録する。負 kinetic 符号の診断は $l=0$ に限定する。

### C.4 Route comparison labels

| label | 内容 | 扱い |
|---|---|---|
| RBD scalar discrepancy | RBD の余分 scalar モード（ $l=0$ で $H_{PP}=-3/8$ ） | 物理的 instability ではない（route-level）。物理主張に不使用 |
| RBD higher-level artifact | $l\ge1$ での RBD extra count | Routh 先行による Lorentz gauge 構造破壊を含むため符号診断は主張しない |
| DBR standard | Dirac 閉包先行の標準 route | 本稿の全主張の基盤 |

### C.5 Scope and limitations (route sector)

- 本稿ではRBD route の出力をポジティブ・ネガティブいずれの物理的主張の根拠としても使用しない。従ってRBD の余分モードを physical ghost / physical instability とは扱わない。
- 本稿では「RBD は常に失敗する」という一般性を主張しない。検証済み sector の対照実験としてのみ扱う。

---

[← Appendix B](DPPUv7-paper07_appB.md) | [Appendix D →](DPPUv7-paper07_appD.md)

# DPPUv7-paper07 — Section 9: Discussion and Outlook

[← Section 8](DPPUv7-paper07_sec08.md) | [Appendix A →](DPPUv7-paper07_appA.md)

---

## 9. Discussion and outlook

### 9.1 What has been established, by layer

本稿の結果は確立の論理的な層に沿って次のように整理できる。以下はすべて DBR プロトコル下の主張である。

**二次レベル**（[Section 4](DPPUv7-paper07_sec04.md), [Section 6](DPPUv7-paper07_sec06.md)）。DBR 二次制約スペクトルにおいて、非伝播 scalar kernel の次元が $d_{\rm zero}=\dim\ker n=3-\mathrm{rank}(n)$ として実現することを証明した。あわせて 2 つの TT モードが正の kinetic を持って伝播すること、 $d_{\rm zero}$ が空間 frame の SO(3) 回転に対して不変であること、および activation rank 機構を示した。本層では kernel を伝播 zero mode と同一視せず、 $d_{\rm zero}$ を新しい不変量とも主張せず、class B への拡張も行わない。

**$k>0$ スペクトル完成**（[Section 5](DPPUv7-paper07_sec05.md)）。解析した harmonic sector において、2 TT の伝播・no-ghost・ $d_{\rm zero}$ への非混入を、 endpoint control（ $T^3$, $S^3$ ）・ $Nil^3$ の閉形式完成・ $Sol^3$ の bulk gap と compact edge の正値性として確認した。これは解析済み harmonic sector に限定した頑健性であり（ $Sol^3$ edge は代表的数値）、full PDE 閉包・all-Anosov 閉形式 edge 定理・compact 連続スペクトルの主張ではない。

**有限変形の分類**（[Section 7](DPPUv7-paper07_sec07.md)）。finite Palatini branch / homogeneous minisuperspace において modulus/wall split（ $T^3$ 3/0, $Nil^3$ 1/1, $Sol^3$ 0/1, $S^3$ 0/0）と minisuperspace no-ghost を示した。これは補助的結果であり、nonlinear lift/mass theorem や moduli 数と $d_{\rm zero}$ の同一視を含まない。

**torsion 閉包**（[Section 8](DPPUv7-paper07_sec08.md)）。検査した AX/VT/MX ansatz の範囲で、reduced torsion sector が exact auxiliary shell に留まり、場昇格でも二次 kinetic を生成せず、有限変形分類が不変であることを示した。これは full 24 成分 $(\partial T)^2$ 不在定理・NY bulk 無関与の一般主張・charge/localization 定理ではない。

これらとは別に、route 比較（[Section 3](DPPUv7-paper07_sec03.md), [Appendix C](DPPUv7-paper07_appC.md)）は対照実験として、RBD route が $S^3$ scalar sector で Dirac 非同値な余分モードを生むことを示す。これは DBR 採用の根拠であり、RBD の一般失敗や physical ghost 解釈を主張するものではない。

### 9.2 新規性の範囲

本稿の topology 階層 $S^3{:}0$ / $Sol^3{:}1$ / $Nil^3{:}2$ / $T^3{:}3$ は、代数的には Bianchi class-A 構造行列の rank-defect・abelianization・ $b_1$ という**既知のデータ**である [16,41,42]。本稿の寄与はその発見ではなく、(i) Lorentzian first-order EC+NY 系で、DBR 制約力学がこのデータを残存する非伝播 scalar kernel として**独立に再現**することの検証、(ii) 単一プロトコル下の four-topology 比較、(iii) その分類が $k>0$ スペクトル完成・有限変形の分類・torsion 閉包の三方向の頑健性検証で保護されること、の三点である。

### 9.3 Relation to Teichmüller/moduli counts of compact Bianchi models

compact Bianchi モデルの global 自由度（Teichmüller パラメータ）は Koike–Tanimoto–Hosoya らによって体系的に勘定されている [24–26]。本稿の $d_{\rm zero}$ はこれと**異なる対象**である: KTH の moduli は metric GR の compact 商の global 構造変形であり、本稿の kernel は first-order EC+NY 系の DBR 制約スペクトルにおける homogeneous 摂動の非伝播方向である。実際、本稿内部でも 有限変形の exact moduli 数（ $T^3$:3, $Nil^3$:1, $Sol^3$:0）は $d_{\rm zero}$（3, 2, 1）と一般に一致しない（[Section 7](DPPUv7-paper07_sec07.md)）。二次 kernel・finite moduli・Teichmüller 変形の三者の精密な対応関係（どの kernel 方向が global 変形に持ち上がるか）は、今後の課題として残る。

### 9.4 Outlook

- **higher-order lift**: finite Palatini branch を超えた full nonlinear lift/mass 解析（wall の力学的な意味の確定）。
- **class-B / VII₀ / VIII への拡張**: 本稿の公式は class A 限定である。class B では対称性 reduction 自体に既知の障害があり [21]、別個の取り扱いを要する。なお Bianchi VII 系のスペクトル問題にも Mathieu 構造が現れることが知られている。
- **$Sol^3$ edge の閉形式化**: representative 数値に留めた compact edge $\lambda_1$ の all-Anosov 閉形式定理への昇格。
- **Euclidean 系列との対応比較**: 本稿の結果（Theorem 1・Table 3・Table 4・Table 5）の、Euclidean response dictionary [4,5] との対応比較。

---

[← Section 8](DPPUv7-paper07_sec08.md) | [Appendix A →](DPPUv7-paper07_appA.md)
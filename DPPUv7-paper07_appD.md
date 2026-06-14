# DPPUv7-paper07 — Appendix D: Endpoint and Spectral Evidence

[← Appendix C](DPPUv7-paper07_appC.md) | [Appendix E →](DPPUv7-paper07_appE.md)

---

## Appendix D. Endpoint and spectral evidence

本 appendix では [Section 4](DPPUv7-paper07_sec04.md)–[5](DPPUv7-paper07_sec05.md) の spectral evidence をまとめる。表記は「解析対象 harmonic sector における $k>0$ スペクトル完成」で統一する。検証基準 H1/H2/H3 は §5 冒頭の operational definition に対応する。

### D.1 Endpoint evidence

- **$T^3$（full- $k$）**: Fourier sector で TT spatial operator $=|k|^2$ を symbolic に確認。 $k=0$ locus は分離（flat degeneracy; scalar-axis zero-mode 全 3 軸は homogeneous レベルの結果）。
- **$S^3$（ $l$-symbolic）**: scalar sector の zero-mode 不在を全 $l$ で symbolic に確認（有限サンプルでない）。spin-2 TT は全 $l\ge2$ で固有値 $(l-1)(l+3)$、kinetic Hessian 正定値。 $l=0,1$ は endpoint 分類（gauge/non-dynamical）。

### D.2 Scalar full-spectral evidence（ $Nil^3$/ $Sol^3$）

非可換 harmonic sector の scalar モードについて、行列調和（matrix harmonic）閉包により有効 zero-mode 数 $d(M)$ を評価:

- $Nil^3$: 非可換モードはすべて正 gap（ $d(M)=0$ ）→ homogeneous 代表の kernel 勘定（ $\dim\ker n=2$ ）が完全な full-spectral 勘定に一致（整合）。
- $Sol^3$: bulk は正 gap（ $d(M)=0$ ）、最下 level の truncation edge 項は定理 zero として数えない（スコープ限定）。

### D.3 TT / spin-1 full-spectral evidence（ $Nil^3$/ $Sol^3$）

同じ matrix harmonic 閉包を spin-2 / spin-1 sector に拡張:

- TT: $d(M)=2$（両 topology）— 伝播 2 TT の full-spectral 整合。
- spin-1: $d(M)=0$（両 topology）— transverse vector に物理自由度なし。

### D.4 Spectral towers（ $k>0$ スペクトル完成の内訳）

| topology | sector | operator / spectrum | grade |
|---|---|---|---|
| $T^3$ | Fourier $k\ne0$ | $\|k\|^2$ | 閉形式（endpoint control） |
| $S^3$ | $l\ge2$ TT | $(l-1)(l+3)$ | 閉形式（endpoint control） |
| $Nil^3$ | toroidal（ $k=0$） | $4\pi^2(m^2+n^2)$ | 閉形式 |
| $Nil^3$ | Heisenberg（ $k\ne0$） | $2\pi\|k\|(2N+1)$ | 閉形式/解析的 |
| $Sol^3$ | boost bulk（ $pq\ne0$） | $\inf V=2\|pq\|>0$ | 解析的 bound（全ラベル） |
| $Sol^3$ | compact edge | 両側 confining modified Mathieu; 離散・ $\lambda_1>0$ [40] | 閉形式/解析的（離散性・正値性）＋ 代表的数値（具体 $\lambda_1$, 代表 $A$ ） |
| $Sol^3$ | 非コンパクト副線（ $pq=0$） | gapless continuum＝amenability の署名 [38,39] | 解釈線（主線と分離; 定理列に不算入） |

### D.5 Spectral verification summary

| 検証基準 | $T^3$ | $S^3$ | $Nil^3$ | $Sol^3$ |
|---|---|---|---|---|
| 伝播 | 閉形式 | 閉形式 | 閉形式/解析的 | 解析的 bound（bulk）+ 代表的数値（edge） |
| no-ghost | 閉形式 | 閉形式 | 閉形式/解析的 | 解析的 bound（bulk; edge TT 符号 positive） |
| 非混入 | 閉形式 | 閉形式 | 閉形式/解析的 | 解析的 bound（bulk; edge carrier 分離） |

### D.6 Independent spectral witness

`script/scripts/paper07/spectral_completion_witness.py` は外部 spectral table を読まず、§5 の式（ $T^3$: $|k|^2$、 $S^3$: $(l-1)(l+3)$、 $Nil^3$ 二系列、 $Sol^3$ bulk gap $\inf V=2|pq|$ と代表 finite-difference edge）をスクリプト内で構成し H1–H3 の正性条件を判定する。

---

[← Appendix C](DPPUv7-paper07_appC.md) | [Appendix E →](DPPUv7-paper07_appE.md)
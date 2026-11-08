# コンパクト斉次幾何上のLorentzian Einstein–Cartan還元：拘束と非斉次摂動

Muacca · **v2** · 改訂日2026-11-08 · 初版2026-06-14

固定コンパクト幾何上のsource-freeなLorentzian Einstein–Cartan理論を、定数完全Nieh–Yan結合のもとで非斉次拘束へ展開する。全対称斉次正則層のT³/Nil/Sol/S³=5/3/2/2対は保持し、対称背景の例外を区別する。flat torus非零modeでは接続保持・計量還元の作用と還元二形式を一致させる。コンパクトCMC初期面などの仮定のもとでArms–Marsden–Moncriefの既知定理を適用し、四つのTaub条件を局所可積分性の必要十分条件とする。単一定在波・非平行二波の明示補正は局所真空方程式を二次まで満たす。時間依存する真空背景上で、LRS Nilの閉じたtensor block、Solの無限dual orbit、S³ Taubの低Wigner blockを構成し、gauge・正準還元・再構成を示す。固定Sol背景ではKID消失時の既知定理により、明示したSobolev級の全実線形化拘束データの局所可積分性を得る。Gowdyと共通coreから再構成したdual-orbitデータを含み、Nil・Taubの二次KID条件は必要条件に留める。有限時間例にはSolの一orbitの独立Lorenz発展と、内部領域での空間・幅・時間の収束対照を含める。4階還元系での収束は未確立であり、全実数直線での数値誤差保証、一般の安定性、固定sector内の非線形閉包、Nil・Taubの持ち上げ十分性は結論しない。

[English README](README.md) / [英語PDF](DPPUv7-paper07.pdf) / [LaTeX](LaTeX/main.tex) / [改訂履歴](CHANGELOG.md) / [参考文献](REFERENCES.md)

## 本文・付録

- [1. 問いと主結果](DPPUv7-paper07_sec01.md)
- [2. 作用・幾何・境界](DPPUv7-paper07_sec02.md)
- [3. 背景と六成分の斉次作用](DPPUv7-paper07_sec03.md)
- [4. 接続消去と拘束保存の対応](DPPUv7-paper07_sec04.md)
- [5. 正則斉次拘束と静的トーラス](DPPUv7-paper07_sec05.md)
- [6. 非斉次摂動の拘束還元](DPPUv7-paper07_sec06.md)
- [7. コンパクト条件と局所可積分性](DPPUv7-paper07_sec07.md)
- [8. 結論](DPPUv7-paper07_sec08.md)
- [付録A. 固定商と規約の辞書](DPPUv7-paper07_appA.md)
- [付録B. 全contortion・Euler式・端点](DPPUv7-paper07_appB.md)
- [付録C. 拘束rankとroute対照の詳細](DPPUv7-paper07_appC.md)
- [付録D. Fourier作用と積分拘束の導出](DPPUv7-paper07_appD.md)
- [付録E. 幾何学的比較資料](DPPUv7-paper07_appE.md)
- [付録F. 真空背景・固定商・調和基底](DPPUv7-paper07_appF.md)
- [付録G. 拘束作用素・還元・再構成](DPPUv7-paper07_appG.md)
- [付録H. 有限時間例とSolの数値対照](DPPUv7-paper07_appH.md)
- [付録I. 改訂注記](DPPUv7-paper07_appI.md)

## 再現用資料

[コード・使用方法](script/README.md)、[保存data](data/README.md)、[採用結果と未確立の範囲](data/perturbations/SCOPE.md)を参照する。Python 3.11以上、repository rootから次を実行する。

```sh
python -m pip install -r script/requirements.txt
python -B script/paper07.py recompute --output results_base_new
python -B script/perturbations.py recompute --model all --output results_models_new
python -B script/paper07.py render --data data/revised --output figures_base_new
python -B script/perturbations.py render --data data/perturbations --output figures_models_new
```

出力先はそれぞれ新しいディレクトリを指定する。model別の実行も可能である。renderは保存dataから描画し、記号導出や数値発展を再実行しない。PDFはLaTeXディレクトリで `pdflatex -interaction=nonstopmode -halt-on-error main.tex` を二回実行する。

斉次正則層、伝播場の種類数、完全な実blockの対数を区別する。flat torusとSolの既知定理の適用、明示二次補正、曲がった背景の線形還元、有限時間数値例は別の根拠を持つ。Solでは整数s≥6の固定商上の初期データを対象に高次の零mode・他orbitを許し、局所時間の存在を結論する。任意のL²還元データや4階還元作用素のwell-posednessは含まない。v1の局所計算engineと保存出力は本補足資料で置き換えた。v1の規約・出力との対応は[比較ガイド](script/legacy_comparison_ja.md)に示す。訂正・撤回の注記は付録IとCHANGELOGにある。

## 配布・謝辞

[Repository](https://github.com/Muacca/DPPUv7-paper07) / [既存Zenodo record](https://doi.org/10.5281/zenodo.20678646)。取得したartifactの版と日付を確認して引用する。[MIT License](LICENSE.txt)。連絡先：[muacca@dmwp.jp](mailto:muacca@dmwp.jp)。

初版ではClaude Opus 4.8、Gemini 3.1 Pro、Grok 4.3、GPT-5.5 Thinking（2026年6月）による支援を謝辞に記した。今回の改訂ではOpenAI Codexが記号導出、コード、日英記述、決定的な図生成を支援した。内容と残る誤りの責任は著者にある。

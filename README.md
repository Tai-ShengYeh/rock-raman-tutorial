# 石膏與大理石的拉曼光譜分析 — R / Python / Orange

> 從化學鍵的振動出發，用三個工具跑完 PCA 與 PLS-DA，再把 VIP 與 loading 一根一根對回化學鍵結。

📖 **[線上閱讀（中文）](https://tai-shengyeh.github.io/rock-raman-tutorial/)** ·
**[Read online (English)](https://tai-shengyeh.github.io/rock-raman-tutorial/index_en.html)**
← 部署後把網址換成你的

*[English summary below ↓](#english)*

---

## 這是什麼

石膏（CaSO₄·2H₂O）與大理石（方解石 CaCO₃）各 5 片、共 10 張拉曼光譜的完整分析教材。
**43 章、68 題互動測驗**，同一套分析用三個工具各做一遍：

| 工具 | 檔案 | 特色 |
|---|---|---|
| **R** | [`r/rock_raman_tutorial.R`](r/rock_raman_tutorial.R)（中文註解）<br>[`r/rock_raman_tutorial_en.R`](r/rock_raman_tutorial_en.R)（English） | 只需 `pls` 一個外掛套件，12 秒跑完 |
| **R（套件版前處理）** | [`r/preprocess_packages.R`](r/preprocess_packages.R)<br>[`r/preprocess_packages_en.R`](r/preprocess_packages_en.R) | `baseline` + `prospectr`，與手寫版逐項比對 |
| **Python** | [`python/rock_raman_tutorial.ipynb`](python/rock_raman_tutorial.ipynb)（中文註解）<br>[`python/rock_raman_tutorial_en.ipynb`](python/rock_raman_tutorial_en.ipynb)（English） | 手寫 NIPALS PLS2，看得見每一步 |
| **Python（套件版前處理）** | [`python/preprocess_packages.ipynb`](python/preprocess_packages.ipynb)<br>[`python/preprocess_packages_en.ipynb`](python/preprocess_packages_en.ipynb) | `pybaselines` + `scipy.signal`，附 arPLS / imodPoly 對照 |
| **Orange** | [`orange/rock_plsda.ows`](orange/rock_plsda.ows) | [orange-spectra](https://pypi.org/project/orange-spectra/) 的 PLS-DA widget，不寫程式 |

**三邊跑出同一組數字**，吻合到 10⁻⁶ 以下 —— 教材第 29 章逐項驗證。

---

## 快速開始

### R
```r
setwd("r")
install.packages("pls")                        # 主分析只需要這一個
source("rock_raman_tutorial.R")                # 約 12 秒

install.packages(c("baseline", "prospectr"))   # 第七篇的套件版前處理
source("preprocess_packages.R")
```

### Python
```bash
pip install numpy scipy pandas matplotlib jupyter
pip install pybaselines                   # 第七篇的套件版前處理
pip install orange-spectra --no-deps      # 選用：第 29 章的對照
jupyter notebook python/rock_raman_tutorial.ipynb
```

### Orange
1. 安裝 [Orange](https://orangedatamining.com/download/)
2. Options → Add-ons → Add more... → `orange-spectra`（或 `pip install orange-spectra`）
3. 開啟 `orange/rock_plsda.ows`，在 File widget 選 `data/rock_preprocessed.tab`
4. PLS-DA widget 的 Components 設 **2**

---

## 檔案結構

```
.
├── index.html                     教材本體・中文（單檔自足，圖已內嵌，離線可讀）
├── index_en.html                  Same tutorial in English
├── data/
│   ├── rock_raman_2class_raw.csv          原始光譜（裁切 150–1800，664 點）
│   ├── rock_preprocessed.csv              扣基線＋平滑後 — 三個工具共用
│   ├── rock_preprocessed.tab              同上，Orange 原生格式（class 已設定）
│   ├── loadings_vip_R.csv                 R 的 PCA loading / PLS loading / VIP
│   ├── qc_table.csv                       各樣品 S/N 與螢光背景
│   ├── jch2018_external_derived.csv       第三方 CC BY 4.0 公開資料（衍生，見 DATA_SOURCES.md）
│   └── rock_raman_3class_raw_including_alabaster.csv   含 alabaster 的完整原始資料
├── r/
│   ├── rock_raman_tutorial.R              中文註解
│   ├── rock_raman_tutorial_en.R           English comments — identical code & output
│   ├── preprocess_packages.R              baseline + prospectr，與手寫版比對
│   └── preprocess_packages_en.R           English
├── python/
│   ├── rock_raman_tutorial.ipynb          中文註解
│   ├── rock_raman_tutorial_en.ipynb       English
│   ├── preprocess_packages.ipynb          pybaselines + scipy，與手寫版比對
│   ├── preprocess_packages_en.ipynb       English
│   └── jch_external_validation.py         外部驗證：JCH 2018 公開資料投影
├── orange/rock_plsda.ows
└── fig/                           17 張輸出圖
```

---

## 主要結果

| | |
|---|---|
| 診斷帶 | 石膏 **1008** cm⁻¹（SO₄ ν₁）／方解石 **1086** cm⁻¹（CO₃ ν₁） |
| PCA（扣基線後） | PC1 = **86.3%**，載荷正向是方解石帶、負向是石膏帶 |
| PLS-DA | 1 個成分即留一 100%，LV1 解釋 86.30% X 變異 |
| 精確排列檢定 | 全部 252 種標籤分派，**p = 2/252 = 0.0079**（＝ 5 對 5 的下限） |
| VIP > 1 | 55 / 664 個波數 |
| ⚠️ VIP 第 2 名 | **419 cm⁻¹ 是儀器假峰** —— 方解石沒有這個振動模式，但兩組都有它 |
| 外部驗證 | 別人的 4 張公開光譜（不同儀器、355 nm 脈衝）投影，校正後 **4/4 判對**，信心 0.82–0.96 |

---

## 教材大綱

1. **化學：為什麼是 1008 和 1086**（7 章）— 極化率、SO₄ 與 CO₃ 的振動模式、為什麼 872 在拉曼看不到、力常數與波數的關係
2. **資料與前處理**（6 章）— 裁切、ALS 基線、Savitzky–Golay、SNV，R 與 Python 並排
3. **R**（5 章）— `prcomp` / `pls::plsr` / VIP / 留一交叉驗證 / 精確排列檢定
4. **Python / Jupyter**（5 章）— SVD、NIPALS 逐行拆解、跨語言重現的坑
5. **Orange**（5 章）— `orange-spectra` PLS-DA widget 的四個輸出、從 VIP 表回到化學鍵
6. **驗證、陷阱與品質把關**（6 章）— 三工具數字對照、儀器假峰、為什麼剔除 alabaster、p ≫ n 的地板、強度響應、18 條錯誤清單
7. **用現成套件做前處理**（4 章）— `baseline` / `prospectr` / `pybaselines` / `scipy.signal`，
   以及三個不會報錯但會安靜改變結果的預設值：`lambda` 的 log10 慣例、SG 砍掉端點、導數的 `delta.wav`
8. **外部驗證**（5 章）— 把 *J. Cult. Herit.* 32 (2018) 的公開資料（355 nm 脈衝 LIF-Raman）
   投影進自己訓練的模型。先抓出資料裡 **46 cm⁻¹ 的波長軸錯誤**（照抄 metadata 的 355 nm，
   PLS-DA 只對 2/4；校正成 354.425 nm 之後 4/4），再談參考文獻、CC BY 4.0 的 attribution 義務，
   以及哪些東西**不能**放上 GitHub

---

## 部署

**GitHub Pages** — Settings → Pages → Source 選 **GitHub Actions**，push 到 `main` 即自動部署
（`.github/workflows/pages.yml` 已備好）。

**Cloudflare Pages** — 連結此 repo，Build command 留空，Build output directory 填 `/`。

兩者都不需要建置步驟，`index.html` 是自足的單一檔案。

---

## 環境

R 4.3.3 · pls 2.8.3 · baseline 1.3-7 · prospectr 0.2.10 · numpy 2.4 · scipy 1.17 · pybaselines 1.2.1 · orange-spectra 0.5.0 · 亂數種子 20260814

## 授權

程式碼 **MIT**；教材文字與圖表 **CC BY-NC-SA 4.0**。

`data/jch2018_external_derived.csv` 是第三方資料的衍生物，原始資料授權 **CC BY 4.0**
（Martínez-Hernández *et al.*, Zenodo, [doi:10.5281/zenodo.3532931](https://doi.org/10.5281/zenodo.3532931)）——
完整出處、授權與修改說明見 **[`DATA_SOURCES.md`](DATA_SOURCES.md)**。

> ⚠️ 相關論文（*J. Cult. Herit.* **32** (2018) 1–8, © Elsevier Masson SAS, All rights reserved）
> **不是**開放取用。論文 PDF 不在本 repo，`.gitignore` 已擋掉 `*.pdf`。
> 論文與其資料集的授權不同 —— 每一次都要各自確認。

帶位歸屬參考礦物拉曼標準光譜，可對照 [RRUFF](https://rruff.info/) 與
[ENS Lyon Handbook of Raman Spectra](https://www.geologie-lyon.fr/Raman/)。


---

<a name="english"></a>

## English

**[Read the tutorial online →](https://tai-shengyeh.github.io/rock-raman-tutorial/index_en.html)**
(`index_en.html` — a single self-contained file; open it locally and it works offline)

Ten Raman spectra — five gypsum (CaSO₄·2H₂O), five marble (calcite, CaCO₃) — analysed three
ways: in **R**, in **Python/Jupyter**, and in **Orange** with the
[orange-spectra](https://pypi.org/project/orange-spectra/) PLS-DA widget.
**43 chapters, 68 interactive questions.**

The point of the tutorial is not the buttons. It is that the same analysis, run in three
completely different interfaces, produces **the same numbers** — and that when the model
returns 100% accuracy it is simultaneously placing second-highest importance on an
**instrument artefact that is not in the sample at all**.

| | |
|---|---|
| Diagnostic bands | gypsum **1008** cm⁻¹ (SO₄ ν₁) / calcite **1086** cm⁻¹ (CO₃ ν₁) |
| PCA (baseline-corrected) | PC1 = **86.3%**; positive loadings are calcite bands, negative are gypsum |
| PLS-DA | leave-one-out 100% with a single component; LV1 explains 86.30% of X variance |
| Exact permutation test | all C(10,5) = 252 label assignments, **p = 2/252 = 0.0079** — the floor for a 5-vs-5 design |
| VIP > 1 | 55 of 664 wavenumbers |
| ⚠️ VIP rank 2 | **419 cm⁻¹ is an instrument artefact** — calcite has no mode there, yet both groups show it |
| ⚠️ LV2 | its four largest loadings *all* fall in that artefact band |

Cross-implementation agreement (Chapter 29): hand-written Python vs the `orange-spectra`
package, **0.000e+00**; Python vs R (`pls::plsr`, `method = "oscorespls"`), **5.0e-07**
(CSV export rounding); the Orange GUI's VIP table matches row by row.

**Contents** — 1. Chemistry: why 1008 and 1086 · 2. Data and preprocessing ·
3. R · 4. Python/Jupyter · 5. Orange · 6. Validation, traps and quality control ·
7. Preprocessing with established packages · 8. External validation

Part 7 swaps the hand-written preprocessing for `baseline` + `prospectr` (R) and
`pybaselines` + `scipy.signal` (Python), and quantifies three defaults that raise no error
but change the result: `lambda` is log10 in R's `baseline` package (a 6% shift in peak
height), `prospectr::savitzkyGolay` discards 4 points at each end (a 12 cm⁻¹ shift in the
VIP table if you forget to trim the axis), and a derivative without `delta.wav` is per
data point rather than per cm⁻¹.

### Quick start

```r
# R
setwd("r"); install.packages("pls"); source("rock_raman_tutorial_en.R")   # ~12 s
```
```bash
# Python
pip install numpy scipy pandas matplotlib jupyter
pip install orange-spectra --no-deps        # optional, for the Chapter 29 comparison
jupyter notebook python/rock_raman_tutorial_en.ipynb
```
```
# Orange
Options → Add-ons → Add more... → orange-spectra
Open orange/rock_plsda.ows, point File at data/rock_preprocessed.tab, Components = 2
```

English-commented versions of both scripts are included:
`r/rock_raman_tutorial_en.R` and `python/rock_raman_tutorial_en.ipynb`.
The code is identical to the Chinese-commented versions and writes byte-for-byte identical
output files.

Part 8 projects four spectra from a published open dataset (Martínez-Hernández *et al.*,
*J. Cult. Herit.* 32, 2018; data on Zenodo under CC BY 4.0) into the model trained here.
Getting there means first catching a 46 cm⁻¹ error in that dataset's wavelength axis: the
metadata's nominal 355 nm displaces every band by 27–58 cm⁻¹, and PLS-DA scores only 2/4.
Refining the excitation wavelength to a single shared 354.425 nm — cross-checked against the
band positions the authors themselves report — gives 4/4 at 0.82–0.96 confidence, and places
the published alabaster with the sulfates, independently confirming the Chapter 31 argument.

Code MIT · text and figures CC BY-NC-SA 4.0 · third-party data CC BY 4.0, see
[`DATA_SOURCES.md`](DATA_SOURCES.md). The related article PDF is all-rights-reserved and is
**not** included.

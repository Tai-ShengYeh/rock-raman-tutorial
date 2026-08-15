# 把 J. Cult. Herit. 32 (2018) 的 alabaster 併進來跑 PCA / PLS-DA

> 結論先講：**可以跑，而且結果對 —— 但要先修一個 46 cm⁻¹ 的波長軸錯誤。**
> 不修的話 PLS-DA 只對 2/4，而且會很有信心地把石灰岩判成石膏。修完 4/4。

## 參考文獻與授權

**論文**
> A. Martínez-Hernández, M. Oujja, M. Sanz, E. Carrasco, V. Detalle, M. Castillejo,
> "Analysis of heritage stones and model wall paintings by pulsed laser excitation of Raman,
> laser-induced fluorescence and laser-induced breakdown spectroscopy signals with a hybrid
> system", *Journal of Cultural Heritage* **32** (2018) 1–8.
> DOI: [10.1016/j.culher.2018.02.004](https://doi.org/10.1016/j.culher.2018.02.004)
> © 2018 Elsevier Masson SAS. **All rights reserved**（非開放取用）

**資料集**
> 同作者群，[Dataset]. Zenodo, 2019-11-08.
> DOI: [10.5281/zenodo.3532931](https://doi.org/10.5281/zenodo.3532931)
> 授權 **CC BY 4.0**．經費：EU IPERION-CH（H2020, Project 654028）

**兩者授權不同。** 資料可以合法重製散布（要標作者、授權、DOI、以及你做了什麼修改），
論文 PDF 不行 —— 見第八節。

**儀器**（論文第 2 節）：Q-switched Nd:YAG（LS-2147, Lotis II）355 nm、17 ns、10 Hz、
6 J cm⁻²；Bentham TMc300 光譜儀，0.30 m，1200 lines/mm 光柵（覆蓋 **50 nm**），狹縫 50 μm；
Andor DH-501 ICCD，delay 0 / gate 50 ns；每張譜累積 125 次。

---

## 一、先看清楚這批檔案是什麼（三個必須先發現的事）

### 1. `Alaba35501–05` 不是五張重複光譜

它們是**同一次量測的五個光柵視窗**，首尾相接：

| 檔案 | 波長範圍 (nm) | 換算 (cm⁻¹) |
|---|---|---|
| Alaba35501 | 362.70 – 415.05 | 598 – 4076 |
| Alaba35502 | 415.96 – 463.85 | 4129 – 6610 |
| Alaba35503 | 463.74 – 517.66 | 6605 – 8851 |
| Alaba35504 | 517.47 – 565.43 | 8844 – 10484 |
| Alaba35505 | 565.16 – 619.78 | 10475 – 12034 |

**只有 `Alaba35501` 涵蓋拉曼指紋區。**其餘四個是 4000–12000 cm⁻¹ 的螢光發射拖尾 ——
xlsx 註明技術是 **LIF-Raman**，這批光譜本來就是螢光為主、拉曼線疊在上面。
所以「五張 alabaster」實際上是 **n = 1**。

### 2. x 軸是奈米，而且 metadata 的「355」是標稱值

`.dat` 的第一欄是波長（nm），不是拉曼位移。要自己換算：

```
Raman shift = 1e7/λ_ex − 1e7/λ
```

照 xlsx 寫的 `355` 去換，四根 ν₁ 帶會落在 963 / 955 / 1028 / 1059 cm⁻¹ ——
**每一根都比該在的位置低了 27–58 cm⁻¹**。實際 λ_ex ≈ **354.425 nm**
（Nd:YAG 三倍頻標稱值就是 354.7 nm，不是 355）。

### 3. `Whitemarble355-001.dat` 與 `-003.dat` 位元組完全相同

MD5 一樣。當成兩個獨立樣品就是自欺欺人。

**能用於指紋區的獨立光譜總數：4 張**（alabaster 1、gypsum 1、limestone 1、white marble 1）。

---

## 二、波長校正 —— 以及它為什麼不是循環論證

用已知的 ν₁ 位置（SO₄ 1008、CO₃ 1086）反推 λ_ex：

| 檔案 | 峰位 (nm) | 指認 | 反推 λ_ex (nm) |
|---|---|---|---|
| Alaba35501 | 367.572 | 1008 | 354.440 |
| YESO35507 | 367.458 | 1008 | 354.333 |
| CALIZA355001 | 368.447 | 1086 | 354.271 |
| Whitemarble355-001 | 368.863 | 1086 | 354.656 |

「用已知峰位反推，再說峰位對了」聽起來像循環論證。兩個檢驗說明它不是：

**檢驗 0（最強的一個）— 論文自己寫了觀測到的峰位。**
第 3.1 節原文：alabaster 與 gypsum 觀測到 197、313、415、494、619、671、**1008**、1135 cm⁻¹
（SO₄²⁻）；limestone 與 marble 觀測到 154、281、711、**1085**、1445、1744 cm⁻¹（CO₃²⁻）。
也就是說 —— **照抄 355.0 換算出來的 963 / 955 / 1028 / 1059，連作者自己報的數字都對不上。**
這是完全獨立於本分析的外部證據，循環論證的疑慮到此為止。

**檢驗 1 — 反推值有沒有按礦物分群？**
硫酸鹽組平均 354.387 nm，碳酸鹽組平均 354.464 nm，差 **77 pm**。
如果指認錯了，兩組會被系統性推開 **979 pm**（= 78 cm⁻¹ 的波長當量）。差了一個數量級。

**檢驗 2 — 對調指認會怎樣？**
把每一根都指給另一個物種，四個 λ_ex 的標準差從 **0.169 nm** 暴增到 **1.187 nm**
（跨度 2.28 nm）—— 等於要假設他們用了四支不同的雷射。

**最關鍵的是**：後面所有分析都只用**一個共用常數** λ_ex = 354.425 nm，不逐檔配。
用這一個常數，四根帶的殘差 RMS = **11.7 cm⁻¹** —— 遠小於這台儀器的
**34–53 cm⁻¹ 頻寬**，也就是說全部落在一個解析度元素之內。

---

## 三、儀器差多少

| | ν₁ 帶 FWHM |
|---|---|
| 你的譜（gypsum 1008） | **12.5 cm⁻¹** |
| 你的譜（marble 1086） | **12.2 cm⁻¹** |
| JCH alabaster | 34 cm⁻¹ |
| JCH gypsum | 48 cm⁻¹ |
| JCH limestone | 43 cm⁻¹ |
| JCH white marble | 53 cm⁻¹ |

**寬了 3–4 倍。** 355 nm 脈衝（17 ns）、125 次累積的 LIF-Raman 混合系統，
色散與狹縫都不是為高解析拉曼設計的。低波數端也整段沒有 —— 公開資料最低到
**598 cm⁻¹**，所以方解石的 156 / 282 晶格模、石膏的 413 / 493 ν₂ 全部拿不到。
共用波數窗只能取 **620–1790 cm⁻¹**。

---

## 四、PCA 結果

19 張譜（你的 15 + 公開的 4），共用軸 620–1790 cm⁻¹ / 5 cm⁻¹ 一點，
前處理完全相同（ALS λ=1e4 → SG 9,3 → SNV）。

| 情境 | PC1 | PC2 | 儀器落差 ÷ 礦物落差（PC1 上） |
|---|---|---|---|
| 照抄 355.0 nm | 37.4% | 26.8% | **0.77** |
| 校正 354.425 nm | 46.5% | 26.9% | **0.17** |
| 校正＋解析度對齊 | 51.0% | 34.5% | **0.18** |

- **照抄**：PC1 上「哪一台儀器」的落差幾乎跟「什麼礦物」一樣大。公開的四張全部擠在
  PC1 負側，不分礦物 —— 分群分的是實驗室。
- **校正後**：PC1 變成乾淨的**硫酸鹽 vs 碳酸鹽**軸，儀器影響掉到礦物影響的 1/6。
- 儀器差異沒有消失，它**搬到 PC2 去了**。公開的四張在 PC2 上一致偏離
  （−5 到 −11），這就是「不同儀器」這件事的座標。**這是好事** —— 混淆因子被隔離在
  一個你看得見、也可以選擇不用的成分上。

把你的譜用 Gaussian 糊到 45 cm⁻¹ FWHM（解析度對齊）幾乎沒有再改善 PC1
（0.17 → 0.18），但把 PC2 拉高到 34.5% —— 說明剩下的儀器差異主要不是頻寬，
而是螢光背景形狀與低波數截斷。

---

## 五、PLS-DA 結果（這才是真正的檢定）

**設計**：只用**你自己的 10 張**（gypsum 5 vs marble 5）訓練，2 個成分；
公開的 4 張完全不參與訓練，當作**外部驗證樣品**投影進去。
這比「全部丟進去一起跑」嚴謹得多 —— 後者在 n=1 的類別上根本沒有意義。

訓練集留一交叉驗證：三種情境都是 **10/10**。

| 樣品 | 真實 | 照抄 355.0 | 校正後 |
|---|---|---|---|
| Alaba35501（alabaster） | 硫酸鹽 | ❌ marble (y_gyp 0.46) | ✅ **gypsum (0.96)** |
| YESO35507（gypsum） | 硫酸鹽 | ✅ gypsum (0.51) | ✅ **gypsum (0.92)** |
| CALIZA355001（limestone） | 碳酸鹽 | ❌ gypsum (0.84) | ✅ **marble (0.82)** |
| Whitemarble355-001 | 碳酸鹽 | ✅ marble (0.46) | ✅ **marble (0.86)** |
| | | **2/4** | **4/4** |

看 y 值比看對錯更有意思：

- **照抄的三個「答案」全部在 0.46–0.51**，也就是模型在丟銅板。唯一有信心的那一個
  （CALIZA 0.84）**是錯的** —— 有信心的錯誤答案，比隨機猜還危險。
- **校正後四個都在 0.82–0.96**，而且方向全對。

### 順帶回答一個舊問題

你自己的 5 張 alabaster，在三種情境下都被判為 **gypsum 5/5**。
公開的那張 alabaster，校正後也是 **gypsum，信心 0.96**。

**兩個實驗室、兩台儀器、兩種激發波長（你的色散式 vs 355 nm 脈衝 LIF-Raman），
獨立得到同一個答案：這兩批 alabaster 都是石膏質的（雪花石膏），不是方解石質的。**

這正是教材第 31 章那個「alabaster 是同一個礦物」論點的外部佐證 —— 而且是別人的資料。

---

## 六、所以能不能用？結論

**能，但只能當外部驗證，不能併進訓練集。** 理由：

1. **n = 1**。alabaster、gypsum、limestone、white marble 各一張獨立光譜。
   PLS-DA 需要類別內變異來估計，n=1 給不出來。硬併進 15 張裡，等於加了 4 個
   高槓桿點，會主導模型而不是驗證模型。
2. **儀器混淆是真的存在的**，只是校正後被推到 PC2。訓練時混進去，模型會學到
   「螢光背景比較高 = JCH = ？」這種與化學無關的規則。
3. **波數窗被截到 620–1790**。丟掉的 156 / 282 / 413 / 493 正是教材裡用來
   佐證礦物指認的低波數證據。

**建議的用法**（也就是本報告做的）：自己的資料訓練，公開資料投影。
這在方法學上是乾淨的，而且比「我的模型 LOO 100%」有力得多 ——
**外部、不同儀器、不同激發波長的樣品也判對了**，這才是泛化能力的證據。

如果你想把它寫進論文或教材，方法段要寫成：

> External validation used four independent spectra from the open dataset of
> Martínez-Hernández *et al.* (*J. Cult. Herit.* 32, 2018), acquired on a hybrid
> LIF–Raman system with 355 nm pulsed excitation (17 ns, 10 Hz, 125 accumulations).
> The published wavelength axis was converted to Raman shift using a single refined
> excitation wavelength of 354.425 nm, obtained from the ν₁ positions of the four
> samples (RMS residual 11.7 cm⁻¹, well within the 34–53 cm⁻¹ instrumental bandwidth);
> taking the nominal 355 nm literally displaces every band by 27–58 cm⁻¹. Spectra were
> interpolated onto a common 620–1790 cm⁻¹ grid at 5 cm⁻¹ and preprocessed identically
> to the training set. All four external samples were correctly classified
> (y ≥ 0.82), including the published alabaster, which was assigned to the
> sulfate class in agreement with our own alabaster specimens.

---

## 七、放上 GitHub：什麼可以、什麼不可以

| 東西 | 授權 | 放上 GitHub |
|---|---|---|
| Zenodo 資料集 `.dat` / `.xlsx` | **CC BY 4.0** | **可以**，須標示出處 |
| Elsevier 論文 PDF `1-s2.0-…-main.pdf` | © 2018 Elsevier Masson SAS，All rights reserved | **不可以** |

最常見的誤會是「補充資料一定跟論文同一個授權」—— 不一定。這篇論文本身在 ScienceDirect 上是
`/article/abs/`（付費牆），版權頁寫 All rights reserved；但作者另外把資料存到 Zenodo 並標 CC BY 4.0。

CC BY 4.0 不是「隨便用」，四件事要做到：**標作者、標授權並連到條款、給原始出處連結、
說明你改了什麼**（我們換算了波數軸、重取樣、只取四個檔）。

已經 commit 過 PDF 的話，光刪檔不夠 —— 歷史裡還在。`git ls-files '*.pdf'` 應該回空。

## 八、還可以做什麼

- `CALIZA355002/003`、`Alaba35502–05`、`Whitemarble355-002` 覆蓋 **2000–12000 cm⁻¹**，
  那是**螢光發射區**。石膏質與方解石質的螢光形狀差很多（見 `jch_overview.png` 上圖），
  可以另外做一組「LIF 而非 Raman」的 PCA —— 這批資料真正獨特的地方其實在這裡。
- `Basemortar532Raman.dat`（灰泥）與 `HgS532Raman.dat`（硃砂）是 **532 nm 的純拉曼**，
  解析度應該好得多，可以當第三個獨立儀器再驗一次。
- `LIBS fresco-266-*.dat` 九張是 LIBS 元素譜，跟拉曼互補（Ca / S / C），
  如果要做「分子＋元素」的資料融合，這批現成可用。

---

## 檔案

- `jch_alabaster_analysis.py` — 完整可重跑的分析（波長校正、反循環檢驗、PCA、PLS-DA 投影）
- `jch_combined.png` — 六格總圖
- `jch_overview.png` — 公開光譜原始樣貌與 900–1200 放大

```bash
python jch_alabaster_analysis.py <JCH 資料夾> rock_raman_3class_raw_including_alabaster.csv
```

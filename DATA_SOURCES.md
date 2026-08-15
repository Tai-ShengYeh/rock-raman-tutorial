# 資料來源與授權 / Data sources and licensing

本 repo 含有三類資料，授權不同，請分開看。
This repository contains three categories of data under different licences.

---

## 1. 自己量測的光譜 / Spectra measured by us

`data/rock_raman_2class_raw.csv`、`data/rock_raman_3class_raw_including_alabaster.csv`、
`data/rock_preprocessed.csv`、`data/rock_preprocessed.tab`、`data/loadings_vip_R.csv`、`data/qc_table.csv`

石膏（gypsum）、大理石（marble）、雪花石膏（alabaster）各 5 片，共 15 張拉曼光譜，
由本 repo 作者量測。**授權：CC BY-NC-SA 4.0**（與教材文字、圖表相同）。

Fifteen Raman spectra (5 gypsum, 5 marble, 5 alabaster) measured by the repository author.
**Licence: CC BY-NC-SA 4.0**, the same as the tutorial text and figures.

---

## 2. 第三方公開資料（衍生） / Third-party open data (derived)

`data/jch2018_external_derived.csv`

**原始出處 / Original source**

> A. Martínez-Hernández, M. Oujja, M. Sanz, E. Carrasco, V. Detalle, M. Castillejo.
> *Analysis of heritage stones and model wall paintings by pulsed laser excitation of Raman,
> laser-induced fluorescence and laser-induced breakdown spectroscopy signals with a hybrid
> system* \[Dataset]. Zenodo, 8 November 2019.
> DOI: [10.5281/zenodo.3532931](https://doi.org/10.5281/zenodo.3532931)

**授權 / Licence**

[Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)

**相關論文 / Related article**

> A. Martínez-Hernández *et al.*, *Journal of Cultural Heritage* **32** (2018) 1–8.
> DOI: [10.1016/j.culher.2018.02.004](https://doi.org/10.1016/j.culher.2018.02.004)
> © 2018 Elsevier Masson SAS. All rights reserved.
>
> ⚠️ **論文 PDF 不在本 repo 中，也不得加入。** 論文與資料集的授權不同：
> 資料是 CC BY 4.0（可重製），論文是全著作權保留。
> **The article PDF is not in this repository and must not be added.** The article and the
> dataset carry different licences: the data is CC BY 4.0, the article is all rights reserved.

**經費 / Funding**

EU IPERION-CH — Integrated Platform for the European Research Infrastructure on Cultural
Heritage (H2020, Project 654028).

### 我們做了什麼修改 / Changes made (required by CC BY 4.0 §3.a.1.B)

1. **波數軸換算。** 原始檔的 x 軸是**波長（nm）**。以 λ_ex = **354.425 nm** 換算成拉曼位移；
   這個值由作者在論文第 3.1 節自己報告的 ν₁ 帶位（SO₄ 1008、CO₃ 1085 cm⁻¹）反推，
   四個檔共用同一個常數，殘差 RMS 11.7 cm⁻¹。
   metadata 標示的 355 nm 是 Nd:YAG 三倍頻的**標稱值**，照抄會讓每根帶低 27–58 cm⁻¹。
   *The x axis was converted from nanometres to Raman shift using a single refined excitation
   wavelength of 354.425 nm, derived from the ν₁ positions the authors report in section 3.1.*
2. **重取樣。** 線性內插到共用格點 620–1790 cm⁻¹，5 cm⁻¹ 一點。
   *Linearly interpolated onto a common 620–1790 cm⁻¹ grid at 5 cm⁻¹.*
3. **只保留 4 個檔。** 原資料集 38 個檔中，只有 4 個涵蓋拉曼指紋區：
   `Alaba35501`、`YESO35507`、`CALIZA355001`、`Whitemarble355-001`。
   `Alaba35502–05` 是**同一次量測的其他光柵視窗**（4000–12000 cm⁻¹ 的 LIF 拖尾）；
   `Whitemarble355-003.dat` 與 `-001.dat` 位元組完全相同，已剔除。
   *Only 4 of 38 files retained. Alaba35502–05 are further grating windows of the same
   measurement; Whitemarble355-003 is byte-identical to -001.*
4. **強度未修改**（任意單位，照原樣）。*Intensities unmodified.*

原始 `.dat` 檔本身**也是 CC BY 4.0，可以合法重製**，但本 repo 選擇只放衍生 CSV 以縮小體積。
要取得原始檔請直接到上面的 Zenodo DOI 下載。
*The original .dat files are also CC BY 4.0 and could legally be redistributed; this repository
ships only the derived CSV to keep it small. Download the originals from the Zenodo DOI above.*

---

## 3. 帶位指認的參考資料 / Reference data for band assignment

教材中的帶位對照（石膏 413/493/620/670/1008/1135、方解石 156/282/712/1086/1435 cm⁻¹）
可對照下列公開資料庫查證。本 repo **未重製**其中任何光譜檔。
The band assignments in the tutorial can be checked against the databases below.
**No spectra from them are reproduced here.**

- RRUFF Project — <https://rruff.info/>
- ENS Lyon, *Handbook of Raman Spectra* — <https://www.geologie-lyon.fr/Raman/>
- IRUG — <http://www.irug.org/>

---

## 引用本 repo / Citing this repository

如果本教材對你的工作有幫助，請同時引用上面第 2 節的原始資料集 ——
CC BY 4.0 的 attribution 義務會隨著資料傳遞下去。
If you use this tutorial, please also cite the original dataset in section 2: the CC BY 4.0
attribution requirement travels with the data.

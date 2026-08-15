# Orange 工作流

## 安裝

1. [下載 Orange](https://orangedatamining.com/download/)
2. 安裝 add-on，兩條路擇一：
   - **Orange 介面**：Options → Add-ons → 右下角 **Add more...** → 輸入 `orange-spectra` → 打勾 → OK → 重開
   - **命令列**：`pip install orange-spectra`（在 Orange 用的那個 Python 環境）

裝好後左側工具列會多一個 **Spectra** 分類，裡面有 11 個 widget。這裡只用 **PLS-DA**。

## 跑起來

1. 開啟 `rock_plsda.ows`
2. 點 **File** widget → 選 `../data/rock_preprocessed.tab`
   - 確認底下 `type` 是 **target**、`ID` 是 **meta**
   - `.tab` 檔已經把這些寫在檔案第三列，正常情況不用改
3. 點 **PLS-DA** → **Components** 設 **2**
4. Status 框應該顯示：
   ```
   10 samples, 2 classes, 2 components.
   Training accuracy: 100.0%
   Confusion (rows = true):
        5    0
        0    5
   ```
5. 點 **VIP table**，第一列應該是 `1086.000`，VIP ≈ 5.718

## 如果 .ows 開不起來

自己拉六個 widget，照這樣連：

```
File ──▶ PLS-DA ──┬──▶ Scatter Plot   (接 Scores)
                  ├──▶ Data Table     (接 VIP)
                  ├──▶ Data Table     (接 Loadings)
                  └──▶ Data Table     (接 Predictions)
```

從 PLS-DA 拉線出去時，Orange 會問你要接哪個輸出 —— 選對應的那個。

## PLS-DA widget 規格

| | |
|---|---|
| 輸入 | **Data**（Table，必須有 discrete class） |
| 輸出 | **Scores**（預設）／**Loadings**／**VIP**（已排序）／**Predictions** |
| 設定 | **Components** 1–20，預設 2 |
| 演算法 | NIPALS PLS2，one-hot 類別，只做 mean-centering、不 scale |
| VIP | `sqrt(p · Σ(ssy_a · (w_ja/‖w_a‖)²) / Σ ssy_a)` |

原始碼：[`orangespectra/widgets/owplsda.py`](https://github.com/Tai-ShengYeh/spectraview/blob/main/orange-spectra/orangespectra/widgets/owplsda.py)
與 [`orangespectra/core.py`](https://github.com/Tai-ShengYeh/spectraview/blob/main/orange-spectra/orangespectra/core.py) 的 `plsda_fit()`

## 注意

Status 框的 **Training accuracy** 是用全部資料建模、再預測同一批資料算出來的。
在 664 個變數對 10 個樣品的情況下它幾乎必然是 100% —— **那個數字不能當成模型效能**。
要交叉驗證請接 Orange 的 **Test and Score**；要精確排列檢定得回到 R 或 Python（教材第 18 / 23 章）。

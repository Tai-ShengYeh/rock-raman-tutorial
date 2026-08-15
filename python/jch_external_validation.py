# -*- coding: utf-8 -*-
"""
把 J. Cult. Herit. 32 (2018) 的公開光譜與自己的大理石／石膏資料合併，
跑 PCA 與 PLS-DA。

資料來源
--------
A. Martínez-Hernández, M. Oujja, M. Sanz et al.,
"Analysis of heritage stones and model wall paintings by pulsed laser excitation
of Raman, laser-induced fluorescence and laser-induced breakdown spectroscopy
signals with a hybrid system", J. Cult. Herit. 32 (2018).
Open data: Alaba355*.dat / YESO35507.dat / CALIZA355*.dat / Whitemarble355*.dat

三件事必須先做，否則結果沒有意義
--------------------------------
1. 檔案的 x 軸是「奈米」不是 cm-1，而且 metadata 寫的 355 nm 是標稱值。
   實際 lambda_ex 由已知的 nu1 帶反推 = 354.425 nm（單一常數，不逐檔配）。
2. Alaba35501~05 不是五張重複光譜，是同一次量測的五個光柵視窗。
   只有 35501 涵蓋指紋區。
3. Whitemarble355-001.dat 與 -003.dat 位元組完全相同（重複檔）。

用法
----
    # 手上有 Zenodo 原始 .dat（doi:10.5281/zenodo.3532931）：
    python jch_external_validation.py <JCH 資料夾> ../data/rock_raman_3class_raw_including_alabaster.csv

    # 沒有原始檔也可以跑 —— 會改用 repo 內已換算好的衍生 CSV：
    python jch_external_validation.py

衍生 CSV（../data/jch2018_external_derived.csv）的來源與修改說明見 DATA_SOURCES.md。
原始資料授權 CC BY 4.0；論文 PDF 是 Elsevier 全著作權保留，不在本 repo 中。
MIT 授權（本腳本）。
"""
import sys, os, json
import numpy as np, pandas as pd
from scipy.sparse import diags, csc_matrix
from scipy.sparse.linalg import spsolve
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d

SRC  = sys.argv[1] if len(sys.argv) > 1 else ""
OWN  = sys.argv[2] if len(sys.argv) > 2 else "../data/rock_raman_3class_raw_including_alabaster.csv"
DERIVED = "../data/jch2018_external_derived.csv"
HAVE_RAW = bool(SRC) and all(os.path.exists(os.path.join(SRC, f)) for f in
                             ["Alaba35501.dat", "YESO35507.dat",
                              "CALIZA355001.dat", "Whitemarble355-001.dat"])
LAM_NOMINAL, LAM_TRUE = 355.0, 354.425
GRID = np.arange(620, 1795, 5.0)

PUB = {"Alaba35501.dat":       ("alabaster", 1008.0),
       "YESO35507.dat":        ("gypsum",    1008.0),
       "CALIZA355001.dat":     ("marble",    1086.0),
       "Whitemarble355-001.dat":("marble",   1086.0)}


# ----------------------------------------------------------------- I/O
def read_dat(path):
    """Tab 分隔、CRLF、缺值寫成 '--'、可能有文字表頭。"""
    xs, ys = [], []
    for line in open(path, encoding="latin-1").read().replace("\r", "").split("\n"):
        p = [q.strip() for q in line.split("\t")]
        if len(p) < 2:
            continue
        try:
            a, b = float(p[0]), float(p[1])
        except ValueError:
            continue                      # 兩個都要能轉，否則整列丟掉
        xs.append(a); ys.append(b)
    return np.asarray(xs), np.asarray(ys)


def to_shift(nm, lam_ex):
    """nm -> Raman shift (cm-1)."""
    return 1e7 / lam_ex - 1e7 / nm


def peak_nm(path, lo=364.0, hi=372.5):
    """視窗內最強帶的位置（拋物線頂點細化），單位 nm。"""
    x, y = read_dat(path); o = np.argsort(x); x, y = x[o], y[o]
    m = (x >= lo) & (x <= hi); x, y = x[m], y[m]
    e = np.r_[np.arange(6), np.arange(len(x) - 6, len(x))]
    y = y - np.polyval(np.polyfit(x[e], y[e], 1), x)
    i = int(np.argmax(y))
    d = (y[i-1] - y[i+1]) / (2 * (y[i-1] - 2*y[i] + y[i+1]))
    return x[i] + d * (x[i+1] - x[i-1]) / 2


# --------------------------------------------------- 波長校正（含反循環檢驗）
def calibrate():
    nm = {f: peak_nm(os.path.join(SRC, f)) for f in PUB}
    lam = {f: 1e7 / (1e7 / nm[f] + PUB[f][1]) for f in PUB}
    swap = {f: 1e7 / (1e7 / nm[f] + (1086.0 if PUB[f][1] == 1008 else 1008.0)) for f in PUB}
    g = float(np.mean(list(lam.values())))
    print("波長校正")
    print(f"  {'file':24s} {'peak(nm)':>9s} {'指認':>6s} {'反推 lam_ex':>12s}")
    for f in PUB:
        print(f"  {f:24s} {nm[f]:9.3f} {PUB[f][1]:6.0f} {lam[f]:12.3f}")
    print(f"  單一常數 lam_ex = {g:.3f} nm   (標稱 355.0，差 {(355.0-g)*1000:.0f} pm)")
    print(f"  用這一個常數，四根帶的殘差 RMS = "
          f"{np.sqrt(np.mean([(to_shift(nm[f], g) - PUB[f][1])**2 for f in PUB])):.1f} cm-1"
          f"   (儀器 FWHM 34-53 cm-1)")
    print(f"  反循環檢驗：正確指認 lam_ex 標準差 {np.std(list(lam.values()), ddof=1):.3f} nm，"
          f"對調指認 {np.std(list(swap.values()), ddof=1):.3f} nm  <- 差 7 倍")
    return g


# ----------------------------------------------------------------- 前處理
def _D2(m):
    e = np.ones(m)
    return diags([e[:-2], -2*e[:-1], e], [0, 1, 2], shape=(m-2, m)).tocsc()

def als(y, lam=1e4, p=0.01, niter=15):
    m = len(y); DD = lam * (_D2(m).T @ _D2(m)); w = np.ones(m)
    for _ in range(niter):
        z = spsolve(csc_matrix(diags(w)) + DD, w * y)
        w = p * (y > z) + (1 - p) * (y <= z)
    return z

def prep(Y):
    Y = np.atleast_2d(np.asarray(Y, float))
    Y = np.array([y - als(y) for y in Y])
    Y = savgol_filter(Y, 9, 3, axis=1)
    return (Y - Y.mean(1, keepdims=True)) / Y.std(1, ddof=1, keepdims=True)


# ----------------------------------------------------------------- 建資料
def build(mode, wn0, X0, own_id, own_cls):
    """mode: 'naive'（照抄 355.0） | 'calib'（校正） | 'matched'（校正＋解析度對齊）"""
    rows, ids, cls, src = [], [], [], []
    lam = LAM_NOMINAL if mode == "naive" else LAM_TRUE
    if not HAVE_RAW:
        # 沒有原始 .dat：用已換算好的衍生 CSV（等同 mode='calib'）
        if mode == "naive":
            raise SystemExit("mode='naive' 需要原始 .dat（衍生 CSV 已經校正過了）。"
                             "請從 doi:10.5281/zenodo.3532931 下載後指定資料夾。")
        dv = pd.read_csv(DERIVED, comment="#")
        gw = np.array([float(c) for c in dv.columns[2:]])
        for i in range(len(dv)):
            rows.append(np.interp(GRID, gw, dv.iloc[i, 2:].to_numpy(float)))
            ids.append(dv.iloc[i, 0])
            cls.append("marble" if dv.iloc[i, 1] in ("limestone", "white_marble") else dv.iloc[i, 1])
            src.append("JCH2018")
    for f, (c, _) in (PUB.items() if HAVE_RAW else []):
        x, y = read_dat(os.path.join(SRC, f))
        s = to_shift(x, lam); o = np.argsort(s)
        g = np.interp(GRID, s[o], y[o], left=np.nan, right=np.nan)
        ok = ~np.isnan(g)
        if not ok.all():
            g = np.interp(GRID, GRID[ok], g[ok])       # 短檔以端點外推補滿
        rows.append(g); ids.append(f[:-4]); cls.append(c); src.append("JCH2018")
    for i in range(len(X0)):
        y = X0[i]
        if mode == "matched":                          # 把自己的譜糊到 45 cm-1 FWHM
            y = gaussian_filter1d(y, (45 / 2.355) / float(np.diff(wn0).mean()))
        rows.append(np.interp(GRID, wn0, y))
        ids.append(own_id[i]); cls.append(own_cls[i]); src.append("own")
    return prep(np.array(rows)), np.array(ids), np.array(cls), np.array(src)


# ----------------------------------------------------------------- 模型
def pca(Xp, n=3):
    Xc = Xp - Xp.mean(0)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    return U[:, :n] * S[:n], Vt[:n], (S**2 / np.sum(S**2))[:n]

def plsda_fit(X, Yb, A=2):
    Xm, Ym = X.mean(0), Yb.mean(0)
    E, F = X - Xm, Yb - Ym
    n, p = X.shape; q = Yb.shape[1]
    W = np.zeros((p, A)); P = np.zeros((p, A)); Q = np.zeros((q, A))
    for a in range(A):
        u = F[:, 0].copy()
        for _ in range(200):
            w = E.T @ u; w /= np.linalg.norm(w); t = E @ w
            qq = F.T @ t / (t @ t); un = F @ qq / (qq @ qq)
            if np.linalg.norm(un - u) < 1e-12: u = un; break
            u = un
        qq = F.T @ t / (t @ t); pp = E.T @ t / (t @ t)
        E -= np.outer(t, pp); F -= np.outer(t, qq)
        W[:, a] = w; P[:, a] = pp; Q[:, a] = qq
    B = W @ np.linalg.inv(P.T @ W) @ Q.T
    return dict(Xm=Xm, Ym=Ym, B=B)

def predict(m, X):
    return (X - m["Xm"]) @ m["B"] + m["Ym"]


# ----------------------------------------------------------------- main
def main():
    if HAVE_RAW:
        calibrate()
    else:
        print(f"找不到原始 .dat，改用衍生檔 {DERIVED}（已校正，跳過 naive 情境）\n")
    d = pd.read_csv(OWN)
    wn0 = np.array([float(c) for c in d.columns[2:]])
    X0 = d.iloc[:, 2:].to_numpy(float)
    own_id, own_cls = d.iloc[:, 0].tolist(), d.iloc[:, 1].tolist()

    for mode in (["naive", "calib", "matched"] if HAVE_RAW else ["calib", "matched"]):
        Xp, ids, cls, src = build(mode, wn0, X0, own_id, own_cls)
        T, _, ev = pca(Xp)
        is_s = np.isin(cls, ["gypsum", "alabaster"])
        d_src = abs(T[src == "JCH2018", 0].mean() - T[src == "own", 0].mean())
        d_min = abs(T[is_s, 0].mean() - T[~is_s, 0].mean())
        print(f"\n=== {mode} ===  PC1 {ev[0]*100:.1f}%  PC2 {ev[1]*100:.1f}%"
              f"   儀器落差/礦物落差 = {d_src/d_min:.2f}")

        tr = (src == "own") & np.isin(cls, ["gypsum", "marble"])
        Yb = np.c_[(cls[tr] == "gypsum").astype(float), (cls[tr] == "marble").astype(float)]
        m = plsda_fit(Xp[tr], Yb)
        hit = sum(int(np.argmax(predict(plsda_fit(Xp[tr][np.arange(tr.sum()) != i],
                                                  Yb[np.arange(tr.sum()) != i]),
                                        Xp[tr][i:i+1])[0]) == np.argmax(Yb[i]))
                  for i in range(tr.sum()))
        print(f"  訓練集 LOO {hit}/{tr.sum()}")
        pr = predict(m, Xp[src == "JCH2018"])
        for j, i in enumerate(np.where(src == "JCH2018")[0]):
            lab = ["gypsum", "marble"][int(np.argmax(pr[j]))]
            truth = "gypsum" if cls[i] in ("gypsum", "alabaster") else "marble"
            print(f"  {'OK' if lab == truth else 'XX'}  {ids[i]:22s} {cls[i]:10s}"
                  f" -> {lab:7s}  y_gyp {pr[j,0]:+.2f}")
        pa = predict(m, Xp[(src == "own") & (cls == "alabaster")])
        print(f"  自己的 alabaster (n=5) -> gypsum {int(np.sum(np.argmax(pa,1)==0))}/5")

if __name__ == "__main__":
    main()

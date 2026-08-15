###############################################################################
#  石膏 vs 大理石　拉曼光譜分析教學腳本（R）
#  PCA + PLS-DA + VIP + loading，全部對回化學鍵結
#
#  資料：data/rock_raman_2class_raw.csv
#        gypsum（石膏 CaSO4.2H2O）5 片 + marble（方解石 CaCO3）5 片
#
#  這份腳本與 python/rock_raman_tutorial.ipynb、以及 Orange 的
#  orange-spectra PLS-DA widget 做完全相同的事，數字吻合到 1e-6 以下。
#
#  只需要一個外掛套件：pls        install.packages("pls")
#  Matrix 是 R 官方隨附，不必另外安裝。
#  程式碼 MIT 授權。
###############################################################################

## ===========================================================================
## 0. 設定
## ===========================================================================
# setwd("path/to/rock-raman-tutorial/r")

if (!requireNamespace("pls", quietly = TRUE)) install.packages("pls")
suppressMessages({ library(Matrix); library(pls) })
set.seed(20260814)
dir.create("../fig", showWarnings = FALSE)

COL <- c(gypsum = "#0F766E", marble = "#4338CA")
PCH <- c(gypsum = 17,        marble = 15)
XL  <- expression(paste("Raman shift (", cm^-1, ")"))


## ===========================================================================
## 1. 讀資料
##    一列一張光譜、一欄一個波數。第 3 欄之後的欄名就是波數。
##    check.names = FALSE 一定要加，否則 "150.004" 會被改成 "X150.004"。
## ===========================================================================
d   <- read.csv("../data/rock_raman_2class_raw.csv", check.names = FALSE)
wn  <- as.numeric(colnames(d)[-(1:2)])
X   <- as.matrix(d[, -(1:2)]);  rownames(X) <- d$ID
grp <- factor(d$type);          ids <- d$ID;  lev <- levels(grp)

cat("光譜張數 =", nrow(X), " 波數點數 =", ncol(X), "\n")
cat("波數範圍 =", round(range(wn), 1), "cm-1\n")
cat("類別 :", paste(lev, table(grp), sep = " x ", collapse = ", "), "\n")
# 664 個波數對 10 個樣品 —— 變數是樣品的 66 倍。這叫 p >> n，第 7 節處理它。


## ===========================================================================
## 2. 前處理三工具
## ===========================================================================

## 2-1 ALS 非對稱最小平方基線：找一條盡量平滑、又盡量從下方貼著資料的線
als <- function(y, lambda = 1e4, p = 0.01, niter = 15) {
  m  <- length(y)
  D  <- diff(Diagonal(m), differences = 2)     # 二階差分 = 曲率
  DD <- lambda * crossprod(D)
  w  <- rep(1, m)
  for (i in 1:niter) {
    W <- Diagonal(x = w)
    z <- as.numeric(solve(W + DD, w * y))
    w <- p * (y > z) + (1 - p) * (y <= z)      # 高於基線者（峰）降權
  }
  z
}

## 2-2 Savitzky-Golay：視窗內配多項式。deriv=1 取一階導數。
sg <- function(y, w = 9, poly = 3, deriv = 0) {
  hw <- (w - 1) / 2; k <- -hw:hw
  A  <- outer(k, 0:poly, "^")
  C  <- solve(crossprod(A), t(A))
  cf <- C[deriv + 1, ] * factorial(deriv)
  n  <- length(y); yp <- c(rep(y[1], hw), y, rep(y[n], hw))
  sapply(1:n, function(i) sum(cf * yp[i:(i + w - 1)]))
}

## 2-3 SNV：每條光譜自己減平均、除標準差，消除整體亮度差異
snv <- function(y) (y - mean(y)) / sd(y)

BASE <- t(apply(X, 1, als))
Xbl  <- t(apply(X - BASE, 1, sg))              # ← 主分析用這個
Xsnv <- t(apply(Xbl, 1, snv))
Xd1  <- t(apply(Xsnv, 1, sg, deriv = 1))
dimnames(Xbl) <- dimnames(Xsnv) <- dimnames(Xd1) <- list(ids, wn)
SETS <- list(Raw = X, Baseline = Xbl, SNV = Xsnv, `SNV+D1` = Xd1)

png("../fig/f02_preprocess.png", 1400, 900, res = 110); par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3, 1))
tt <- c("A. Raw", "B. + ALS baseline + SG smooth", "C. + SNV", "D. + 1st derivative")
for (k in 1:4) {
  matplot(wn, t(SETS[[k]]), type = "l", lty = 1, lwd = 1.3, col = COL[as.character(grp)],
          xlab = XL, ylab = "", main = tt[k])
  abline(v = c(1008, 1086), col = "grey70", lty = 3)
  if (k == 1) legend("topright", lev, col = COL, lwd = 2, bty = "n", cex = .85)
}
dev.off()


## ===========================================================================
## 3. 峰位鑑定：先做化學，再做統計
## ---------------------------------------------------------------------------
##  石膏  CaSO4.2H2O : 413/493 (SO4 v2)、620/670 (v4)、1008 (v1 最強)、1135 (v3)
##  方解石 CaCO3      : 156/282 (晶格)、712 (CO3 v4)、1086 (v1 最強)、1435 (v3)
## ===========================================================================
BANDS <- data.frame(
  wn  = c(156, 282, 413, 493, 620, 670, 712, 1008, 1086, 1135, 1435),
  who = c("cal","cal","gyp","gyp","gyp","gyp","cal","gyp","cal","gyp","cal"))

peakh <- function(M, t, tol = 8) {
  j <- which(wn >= t - tol & wn <= t + tol); apply(M[, j, drop = FALSE], 1, max) }
snr <- function(M, t) {
  bg <- M[, (wn > 930 & wn < 985) | (wn > 1150 & wn < 1250), drop = FALSE]
  (peakh(M, t) - rowMeans(bg)) / apply(bg, 1, sd) }

QC <- data.frame(sample = ids, type = as.character(grp),
                 SN_1008 = round(snr(Xbl, 1008), 1),
                 SN_1086 = round(snr(Xbl, 1086), 1),
                 h419    = round(peakh(Xbl, 419, 4), 1),
                 fluor   = round(rowMeans(BASE), 1))
print(QC, row.names = FALSE)
write.csv(QC, "../data/qc_table.csv", row.names = FALSE)
# 門檻：S/N >= 3 才算偵測到，>= 10 才算可靠。兩組的診斷帶都遠超門檻。
# 但注意 h419 那一欄 —— 419 在兩組都有，而方解石沒有這個振動模式。儀器假峰。


## ===========================================================================
## 4. PCA（非監督式）
##    center = TRUE 一定要；scale. = FALSE —— 光譜不要標準化每個波數欄位，
##    否則只有雜訊的波數會被放大到跟 1086 主帶一樣重要。
## ===========================================================================
cat("\nPCA 解釋變異 (%)：\n")
for (nm in names(SETS)) {
  p <- prcomp(SETS[[nm]], center = TRUE, scale. = FALSE)
  v <- 100 * p$sdev^2 / sum(p$sdev^2)
  cat(sprintf("  %-9s PC1 %5.1f  PC2 %5.1f  PC3 %5.1f\n", nm, v[1], v[2], v[3]))
}
pc <- prcomp(Xbl, center = TRUE, scale. = FALSE)
ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)

png("../fig/f06_pca_scores.png", 1400, 900, res = 110); par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3, 1))
for (nm in names(SETS)) {
  p <- prcomp(SETS[[nm]], center = TRUE, scale. = FALSE)
  v <- 100 * p$sdev^2 / sum(p$sdev^2); s <- p$x
  plot(s[, 1], s[, 2], col = COL[as.character(grp)], pch = PCH[as.character(grp)], cex = 1.9,
       xlab = sprintf("PC1 (%.1f%%)", v[1]), ylab = sprintf("PC2 (%.1f%%)", v[2]),
       main = paste("PCA -", nm))
  abline(h = 0, v = 0, col = "grey85", lty = 3)
  for (g in lev) { k <- grp == g; h <- chull(s[k, 1:2])       # 凸包，n=5 不畫信賴橢圓
    polygon(s[k, 1][c(h, h[1])], s[k, 2][c(h, h[1])],
            border = COL[g], col = adjustcolor(COL[g], .10), lty = 2) }
  if (nm == "Raw") legend("topright", lev, col = COL, pch = PCH, bty = "n", cex = .85)
}
dev.off()

## 載荷圖：得分圖說「誰跟誰像」，載荷圖才說「因為哪根峰」
cat("\nPC1 |loading| 最大的 10 個波數：\n")
L <- pc$rotation[, 1]; o <- order(-abs(L))[1:10]
print(data.frame(wn = round(wn[o], 1), loading = round(L[o], 4),
                 side = ifelse(L[o] > 0, "calcite", "gypsum")), row.names = FALSE)


## ===========================================================================
## 5. PLS-DA（監督式）
##    把類別編成 0/1 假變數矩陣，對 X 做 PLS2 迴歸，取最大輸出當預測類別。
##    method = "oscorespls" 就是 NIPALS PLS2 —— 與 orange-spectra 的
##    plsda_fit() 同一套演算法；simpls 不回傳 loading.weights，算不了 VIP。
## ===========================================================================
Ydum    <- function(g) { m <- model.matrix(~ g - 1); colnames(m) <- levels(g); m }
fitpls  <- function(X, g, nc) plsr(Ydum(g) ~ X, ncomp = nc, method = "oscorespls", scale = FALSE)
predcls <- function(m, Xn, lev, nc) {
  pr <- predict(m, newdata = Xn, ncomp = nc)
  factor(lev[apply(matrix(pr, ncol = length(lev)), 1, which.max)], levels = lev) }

mod <- suppressWarnings(fitpls(Xbl, grp, 4))
cat("\n各成分解釋的 X 變異 (%)：\n"); print(round(pls::explvar(mod), 2))

png("../fig/f09_plsda_scores.png", 820, 640, res = 110); par(mar = c(4.2, 4.2, 3, 1))
s <- mod$scores[, 1:2]; ev <- pls::explvar(mod)
plot(s, col = COL[as.character(grp)], pch = PCH[as.character(grp)], cex = 2,
     xlab = sprintf("LV1 (X-var %.2f%%)", ev[1]), ylab = sprintf("LV2 (X-var %.2f%%)", ev[2]),
     main = "PLS-DA scores")
abline(h = 0, v = 0, col = "grey85", lty = 3)
text(s[, 1], s[, 2], sub("[a-z]+", "", ids), cex = .65, pos = 1)
legend("topright", lev, col = COL, pch = PCH, bty = "n"); dev.off()


## ===========================================================================
## 6. VIP：模型在看哪些波數
##    VIP_j = sqrt( p * sum_a( ssy_a * (w_ja/||w_a||)^2 ) / sum_a ssy_a )
##    VIP > 1 是慣用門檻（VIP 平方的平均剛好等於 1）。
##    這條公式與 orange-spectra core.py 裡的完全相同。
## ===========================================================================
vip <- function(m, nc) {
  W  <- m$loading.weights[, 1:nc, drop = FALSE]
  T  <- m$scores[,          1:nc, drop = FALSE]
  Q  <- m$Yloadings[,       1:nc, drop = FALSE]
  SS <- colSums(Q^2) * colSums(T^2)              # 每個成分解釋掉的 Y 變異
  Wn <- sweep(W, 2, sqrt(colSums(W^2)), "/")
  sqrt(nrow(W) * rowSums(sweep(Wn^2, 2, SS, "*")) / sum(SS))
}
NC <- 2
V  <- vip(mod, NC)
cat("\nVIP > 1 的波數：", sum(V > 1), "/", length(V), "\n")
cat("VIP 前 15 名：\n")
print(head(data.frame(wn = round(wn, 1), VIP = round(V, 3))[order(-V), ], 15), row.names = FALSE)
# 第 2、5、7 名的 419 / 422 / 416 是儀器假峰 —— VIP 高只代表「模型用了它」。

png("../fig/f10_vip.png", 1500, 860, res = 110); par(mfrow = c(2, 1), mar = c(0.4, 4.6, 3, 1))
GM <- t(sapply(lev, function(g) colMeans(Xbl[grp == g, , drop = FALSE])))
matplot(wn, t(GM), type = "l", lty = 1, lwd = 1.9, col = COL, xaxt = "n", ylab = "Intensity",
        main = "VIP > 3 vs group mean spectra")
abline(v = wn[V > 3], col = adjustcolor("#F59E0B", .30), lwd = 2)
matlines(wn, t(GM), lty = 1, lwd = 1.9, col = COL)
legend("topright", lev, col = COL, lwd = 2, bty = "n", cex = .85)
par(mar = c(4.2, 4.6, 0.4, 1))
plot(wn, V, type = "h", col = ifelse(V > 3, "#B45309", ifelse(V > 1, "grey55", "grey82")),
     xlab = XL, ylab = "VIP")
abline(h = 1, col = "#B91C1C", lty = 2)
abline(v = BANDS$wn, col = "grey88", lty = 3)
text(BANDS$wn, max(V) * .93, BANDS$wn, srt = 90, cex = .6,
     col = ifelse(BANDS$who == "gyp", COL["gypsum"], COL["marble"]))
dev.off()


## ===========================================================================
## 7. 驗證
## ===========================================================================

## 7-1 留一交叉驗證：每一輪都要連模型一起重建
loocv <- function(M, g, ncmax = 4) {
  n <- nrow(M); res <- matrix(NA, n, ncmax)
  for (i in 1:n) {
    m <- fitpls(M[-i, , drop = FALSE], factor(g[-i], levels = lev), ncmax)
    for (a in 1:ncmax) res[i, a] <- as.character(predcls(m, M[i, , drop = FALSE], lev, a))
  }
  apply(res, 2, function(p) mean(p == as.character(g)))
}
cat("\nPLS-DA 留一正確率：\n")
for (nm in names(SETS))
  cat(sprintf("  %-9s %s\n", nm,
      paste(sprintf("%dLV=%.3f", 1:4, suppressWarnings(loocv(SETS[[nm]], grp))), collapse = "  ")))

pred <- suppressWarnings(sapply(1:10, function(i) {
  m <- fitpls(Xbl[-i, , drop = FALSE], factor(grp[-i], levels = lev), NC)
  as.character(predcls(m, Xbl[i, , drop = FALSE], lev, NC)) }))
cm <- table(true = grp, predicted = factor(pred, levels = lev))
cat("\n留一混淆矩陣（2 LV）：\n"); print(cm)
obs <- mean(pred == as.character(grp)); cat("正確率 =", obs, "\n")

## 7-2 精確排列檢定
##  5 對 5 只有 C(10,5) = 252 種標籤分派 —— 全部跑完就得到精確 p 值，不必抽樣。
cmb  <- combn(10, 5)
accs <- suppressWarnings(apply(cmb, 2, function(idx) {
  gg <- factor(ifelse(1:10 %in% idx, lev[1], lev[2]), levels = lev)
  mean(sapply(1:10, function(i) {
    m <- fitpls(Xbl[-i, , drop = FALSE], factor(gg[-i], levels = lev), NC)
    as.character(predcls(m, Xbl[i, , drop = FALSE], lev, NC)) }) == as.character(gg)) }))
pval <- sum(accs >= obs) / length(accs)
cat(sprintf("\n精確排列檢定：觀測 = %.3f，隨機平均 = %.3f，p = %d/%d = %.4f\n",
            obs, mean(accs), sum(accs >= obs), length(accs), pval))
# 只有 2 種分派達到 100%：真實標籤，以及把兩組完全對調的那一種。
# 所以 p = 2/252 = 0.0079 是這個樣本數的下限，不可能更小。

png("../fig/f11_perm.png", 900, 600, res = 110); par(mar = c(4.2, 4.2, 3, 1))
h <- hist(accs, breaks = seq(-0.05, 1.05, by = 0.1), col = "#CBD5E1", border = "white",
          xlab = "LOO-CV accuracy under permuted labels",
          main = sprintf("Exact permutation test: all %d label assignments", length(accs)))
abline(v = obs, col = "#B91C1C", lwd = 3)
text(obs, max(h$counts) * .75, sprintf("observed = %.2f\np = %d/%d = %.4f ",
     obs, sum(accs >= obs), length(accs), pval), col = "#B91C1C", adj = 1, cex = .85)
dev.off()


## ===========================================================================
## 8. 匯出：給 Python / Orange 對照用
## ===========================================================================
out <- data.frame(ID = ids, type = as.character(grp), round(Xbl, 6), check.names = FALSE)
colnames(out) <- c("ID", "type", format(wn, trim = TRUE))
write.csv(out, "../data/rock_preprocessed.csv", row.names = FALSE)

write.csv(data.frame(wn = round(wn, 4), VIP = round(V, 6),
                     PC1 = round(pc$rotation[, 1], 6), PC2 = round(pc$rotation[, 2], 6),
                     PLS_p1 = round(mod$loadings[, 1], 6),
                     PLS_p2 = round(mod$loadings[, 2], 6)),
          "../data/loadings_vip_R.csv", row.names = FALSE)

## Orange 原生 .tab（第 1 列名稱、第 2 列型別、第 3 列旗標）
tab <- "../data/rock_preprocessed.tab"
cat(paste(c("ID", format(wn, trim = TRUE), "type"), collapse = "\t"), "\n", sep = "", file = tab)
cat(paste(c("string", rep("c", length(wn)), "d"), collapse = "\t"), "\n", sep = "", file = tab, append = TRUE)
cat(paste(c("meta", rep("", length(wn)), "class"), collapse = "\t"), "\n", sep = "", file = tab, append = TRUE)
for (i in 1:nrow(Xbl))
  cat(paste(c(ids[i], sprintf("%.6f", Xbl[i, ]), as.character(grp[i])), collapse = "\t"),
      "\n", sep = "", file = tab, append = TRUE)

cat("\n完成。圖在 ../fig/，資料在 ../data/。\n")
###############################################################################

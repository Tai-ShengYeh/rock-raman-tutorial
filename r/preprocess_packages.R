###############################################################################
#  用現成套件做前處理：baseline 與 prospectr
#  ---------------------------------------------------------------------------
#  主教材 rock_raman_tutorial.R 裡的 ALS / Savitzky-Golay / SNV 都是手寫的，
#  目的是讓你看見每一步在算什麼。實務上你會用套件。這份腳本把兩邊擺在一起跑，
#  然後逐項比對數字 —— 並且把三個「預設值會咬人」的地方標出來。
#
#      install.packages(c("baseline", "prospectr"))
#
#  這兩個套件都不是必要的：沒裝的話，腳本會改用「照 CRAN 原始碼逐行複刻」的
#  版本繼續跑完，並在輸出中標明。裝了才是真正的驗證。
#
#  資料：../data/rock_raman_2class_raw.csv
#  MIT 授權。
###############################################################################

suppressMessages(library(Matrix))
HAS_BASELINE  <- requireNamespace("baseline",  quietly = TRUE)
HAS_PROSPECTR <- requireNamespace("prospectr", quietly = TRUE)
cat("baseline  :", if (HAS_BASELINE)  "installed" else "NOT installed (using replication)", "\n")
cat("prospectr :", if (HAS_PROSPECTR) "installed" else "NOT installed (using replication)", "\n\n")

d   <- read.csv("../data/rock_raman_2class_raw.csv", check.names = FALSE)
wn  <- as.numeric(colnames(d)[-(1:2)])
X   <- as.matrix(d[, -(1:2)]);  rownames(X) <- d$ID
STEP <- mean(diff(wn))                       # 2.4863 cm-1
cat(sprintf("X = %d x %d,  wn %.1f-%.1f,  mean step %.4f cm-1\n\n",
            nrow(X), ncol(X), min(wn), max(wn), STEP))


## ===========================================================================
## 1. 手寫版（與主教材完全相同，抄過來當基準）
## ===========================================================================
als_hand <- function(y, lambda = 1e4, p = 0.01, niter = 15) {
  m <- length(y); D <- diff(Diagonal(m), differences = 2); DD <- lambda * crossprod(D)
  w <- rep(1, m)
  for (i in 1:niter) {
    z <- as.numeric(solve(Diagonal(x = w) + DD, w * y))
    w <- p * (y > z) + (1 - p) * (y <= z)
  }
  z
}
sg_hand <- function(y, w = 9, poly = 3, deriv = 0) {
  hw <- (w - 1) / 2; k <- -hw:hw
  A  <- outer(k, 0:poly, "^"); C <- solve(crossprod(A), t(A))
  cf <- C[deriv + 1, ] * factorial(deriv)
  n  <- length(y); yp <- c(rep(y[1], hw), y, rep(y[n], hw))   # 邊界：複製端點
  sapply(1:n, function(i) sum(cf * yp[i:(i + w - 1)]))
}
snv_hand <- function(y) (y - mean(y)) / sd(y)


## ===========================================================================
## 2. 複刻版：照 CRAN 原始碼一行一行搬過來，沒裝套件時當備援
## ---------------------------------------------------------------------------
##  這兩段不是「我覺得套件應該這樣做」，是直接對照
##    https://github.com/cran/baseline/blob/master/R/baseline.als.R
##    https://github.com/cran/prospectr/blob/master/R/savitzkyGolay.R
##  寫出來的。註解裡標出與手寫版的每一個差異。
## ===========================================================================

## baseline::baseline.als 的核心
##   差異 1：lambda 是「10 的次方」—— DD <- .create_DD(m) * 10^lambda
##   差異 2：權重更新用 (y < z) 而非 (y <= z)，等號落在 z 上的點權重是 0
##   差異 3：權重不再變動就提早跳出（手寫版是固定跑滿 niter）
##   差異 4：預設 p = 0.05（手寫版用 0.01）
als_repl <- function(y, lambda = 6, p = 0.05, maxit = 20) {
  m <- length(y); D <- diff(Diagonal(m), differences = 2)
  DD <- 10^lambda * crossprod(D)
  w <- rep(1, m)
  for (it in 1:maxit) {
    z <- as.numeric(solve(Diagonal(x = w) + DD, w * y))
    w_old <- w
    w <- p * (y > z) + (1 - p) * (y < z)
    if (sum(w_old != w) == 0) break
  }
  z
}

## prospectr::savitzkyGolay 的核心
##   差異 1：valid convolution —— 輸出比輸入短 (w-1) 點，兩端各砍 (w-1)/2
##   差異 2：有 delta.wav 參數，導數會除以 delta.wav^m 變成「每 cm-1」
sg_repl <- function(y, m = 0, p = 3, w = 9, delta.wav = NULL) {
  gap   <- (w - 1) / 2
  basis <- outer(-gap:gap, 0:p, "^")
  A     <- solve(crossprod(basis, basis), t(basis))
  cf    <- factorial(m) * A[m + 1, ]
  n <- length(y)
  out <- sapply((gap + 1):(n - gap), function(i) sum(cf * y[(i - gap):(i + gap)]))
  if (!is.null(delta.wav)) out <- out / delta.wav^m
  out
}


## ===========================================================================
## 3. 真正呼叫套件（有裝才跑）
## ===========================================================================

## ---- baseline ----------------------------------------------------------
## 注意三件事：
##   (a) spectra 必須是 matrix，單張光譜要 drop = FALSE，否則 baseline() 會報錯
##   (b) 光譜放在「列」（rows），跟本教材的排法一致
##   (c) lambda = 4 才等於手寫版的 lambda = 1e4
als_pkg <- function(Y, lambda = 4, p = 0.01, maxit = 20) {
  if (!is.matrix(Y)) Y <- matrix(Y, nrow = 1)
  bc <- baseline::baseline(Y, method = "als", lambda = lambda, p = p, maxit = maxit)
  baseline::getBaseline(bc)          # getCorrected(bc) 直接給扣完的譜
}

## ---- prospectr ---------------------------------------------------------
## savitzkyGolay(X, m, p, w, delta.wav)
##   m = 導數階數, p = 多項式階數, w = 視窗長度（奇數）
## standardNormalVariate(X) 逐列做，分母是 n-1（跟 R 的 sd() 一樣）
sg_pkg  <- function(Y, m = 0, p = 3, w = 9, delta.wav = NULL) {
  if (is.null(delta.wav)) prospectr::savitzkyGolay(Y, m = m, p = p, w = w)
  else prospectr::savitzkyGolay(Y, m = m, p = p, w = w, delta.wav = delta.wav)
}
snv_pkg <- function(Y) prospectr::standardNormalVariate(Y)


## ===========================================================================
## 4. 逐項比對
## ===========================================================================
line <- function(ch = "-") cat(strrep(ch, 74), "\n")
rmsd <- function(a, b) sqrt(mean((a - b)^2))
ok   <- function(v, tol = 1e-8) if (v < tol) "  <= identical" else ""

y  <- X[1, ]
gap <- 4                                    # (w-1)/2，w = 9

line("="); cat("A. ALS 基線\n"); line()

zh  <- als_hand(y)                          # lambda = 1e4
zr4 <- als_repl(y, lambda = 4, p = 0.01)    # 複刻，lambda = 4 -> 10^4
zr6 <- als_repl(y)                          # 複刻，套件預設 lambda = 6, p = .05
cat(sprintf("  hand(1e4)  vs  replication(lambda=4, p=.01) : RMSD %.3e%s\n",
            rmsd(zh, zr4), ok(rmsd(zh, zr4))))
cat(sprintf("  hand(1e4)  vs  replication(lambda=6, p=.05) : RMSD %8.4f counts\n",
            rmsd(zh, zr6)))
if (HAS_BASELINE) {
  zp4 <- as.numeric(als_pkg(y, lambda = 4, p = 0.01))
  zp6 <- as.numeric(als_pkg(y, lambda = 6, p = 0.05))
  cat(sprintf("  hand(1e4)  vs  baseline::als(lambda=4,p=.01): RMSD %.3e%s   <-- REAL PACKAGE\n",
              rmsd(zh, zp4), ok(rmsd(zh, zp4))))
  cat(sprintf("  replication vs baseline::als(lambda=6,p=.05): RMSD %.3e%s   <-- REAL PACKAGE\n",
              rmsd(zr6, zp6), ok(rmsd(zr6, zp6))))
}
pk <- function(z, c = 1008, tol = 8) max((y - z)[wn > c - tol & wn < c + tol])
cat(sprintf("\n  1008 cm-1 峰高： hand %.1f | lambda=4 %.1f | lambda=6(預設) %.1f  -> 差 %.1f%%\n",
            pk(zh), pk(zr4), pk(zr6), 100 * (pk(zh) - pk(zr6)) / pk(zh)))
cat("  >> lambda 是 log10。照抄手寫版的 1e4 填進去會變成 10^10000。\n")

line("="); cat("B. Savitzky-Golay\n"); line()

h  <- sg_hand(y)
r  <- sg_repl(y)
cat(sprintf("  長度： hand %d | prospectr 複刻 %d | 輸入 %d   (少了 w-1 = 8 點)\n",
            length(h), length(r), length(y)))
cat(sprintf("  內部（兩端各去 %d 點） hand vs 複刻 : max diff %.3e%s\n",
            gap, max(abs(h[(gap + 1):(length(y) - gap)] - r)),
            ok(max(abs(h[(gap + 1):(length(y) - gap)] - r)))))
if (HAS_PROSPECTR) {
  pp <- as.numeric(sg_pkg(y))
  cat(sprintf("  內部  hand vs prospectr::savitzkyGolay : max diff %.3e%s   <-- REAL PACKAGE\n",
              max(abs(h[(gap + 1):(length(y) - gap)] - pp)),
              ok(max(abs(h[(gap + 1):(length(y) - gap)] - pp)))))
}
cat(sprintf("  但邊界不一樣： hand 的第 1 點 = %.2f，prospectr 根本沒有這一點。\n", h[1]))
cat("  >> 用 prospectr 之後波數軸要跟著砍：wn_sg <- wn[5:(length(wn)-4)]\n")

cat(sprintf("\n  一階導數的尺度： 不給 delta.wav  max|d1| = %.4f  (每「點」)\n",
            max(abs(sg_repl(y, m = 1)))))
cat(sprintf("                   給 delta.wav   max|d1| = %.4f  (每 cm-1，差 %.4f 倍)\n",
            max(abs(sg_repl(y, m = 1, delta.wav = STEP))), STEP))
cat("  >> 只做分類不影響（全體同乘一個常數），要報導數數值就一定要給。\n")

line("="); cat("C. SNV\n"); line()

A <- t(apply(X, 1, snv_hand))
B <- { C0 <- X - rowMeans(X); C0 / sqrt(rowSums(C0^2) / (ncol(X) - 1)) }
cat(sprintf("  hand vs prospectr 複刻 : max diff %.3e%s\n", max(abs(A - B)), ok(max(abs(A - B)))))
if (HAS_PROSPECTR) {
  P <- snv_pkg(X)
  cat(sprintf("  hand vs standardNormalVariate : max diff %.3e%s   <-- REAL PACKAGE\n",
              max(abs(A - P)), ok(max(abs(A - P)))))
}
Cn <- (X - rowMeans(X)) / apply(X, 1, function(r) sqrt(mean((r - mean(r))^2)))
cat(sprintf("  如果分母用 n 而不是 n-1 : max diff %.3e  (差 sqrt(n/(n-1)) = %.6f 倍)\n",
            max(abs(A - Cn)), sqrt(ncol(X) / (ncol(X) - 1))))
cat("  >> 這一個手寫與套件完全一致，因為兩邊都用 n-1。\n")


## ===========================================================================
## 5. 實務寫法：整條前處理鏈，套件版
## ===========================================================================
line("="); cat("D. 套件版的完整前處理鏈\n"); line()

if (HAS_BASELINE && HAS_PROSPECTR) {
  bc    <- baseline::baseline(X, method = "als", lambda = 4, p = 0.01, maxit = 20)
  Xcorr <- baseline::getCorrected(bc)                       # 10 x 664
  Xsg   <- prospectr::savitzkyGolay(Xcorr, m = 0, p = 3, w = 9)   # 10 x 656
  wn_sg <- wn[5:(length(wn) - 4)]                           # 波數軸跟著砍
  Xsnv  <- prospectr::standardNormalVariate(Xsg)
  Xd1   <- prospectr::savitzkyGolay(Xsnv, m = 1, p = 3, w = 9, delta.wav = STEP)
  wn_d1 <- wn_sg[5:(length(wn_sg) - 4)]                     # 再砍一次！
  cat(sprintf("  Xcorr %d x %d -> Xsg %d x %d -> Xd1 %d x %d\n",
              nrow(Xcorr), ncol(Xcorr), nrow(Xsg), ncol(Xsg), nrow(Xd1), ncol(Xd1)))
  cat(sprintf("  波數軸： %d -> %d -> %d\n", length(wn), length(wn_sg), length(wn_d1)))
  cat("  >> 每過一次 savitzkyGolay 就少 8 個變數。做兩次導數鏈就少 16 個。\n")
  cat("  >> VIP 表要對回波數時，用的必須是砍過的 wn_d1，不是原始的 wn。\n")
  write.csv(data.frame(wavenumber = wn_sg, t(Xsg)),
            "../data/rock_preprocessed_pkg.csv", row.names = FALSE)
  cat("  已寫出 ../data/rock_preprocessed_pkg.csv\n")
} else {
  cat("  兩個套件都裝了才會跑這一段。\n")
  cat("  install.packages(c(\"baseline\", \"prospectr\"))\n")
}

line("="); cat("結論\n"); line()
cat("  1. lambda 在 baseline 套件裡是 log10 —— 手寫的 1e4 要寫成 lambda = 4。\n")
cat("  2. prospectr::savitzkyGolay 會砍掉兩端各 (w-1)/2 點，波數軸要同步砍。\n")
cat("  3. SNV 兩邊一致（都是 n-1），這一個可以放心換。\n")
cat("  4. 對回化學鍵之前，先確認你的波數向量跟資料矩陣還是對齊的。\n")

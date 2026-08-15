###############################################################################
#  Preprocessing with established packages: baseline and prospectr
#  ---------------------------------------------------------------------------
#  In the main script rock_raman_tutorial.R the ALS / Savitzky-Golay / SNV
#  routines are hand-written, so that you can see what every step computes.
#  In practice you would use a package. This script runs both side by side and
#  compares them item by item -- and marks the three defaults that will bite you.
#
#      install.packages(c("baseline", "prospectr"))
#
#  Neither package is required: without them the script falls back to a
#  line-by-line replication of the CRAN sources and says so in the output.
#  Installing them is what turns this into a real verification.
#
#  Data: ../data/rock_raman_2class_raw.csv
#  MIT licence.
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
## 1. Hand-written (identical to the main script; the reference here)
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
  n  <- length(y); yp <- c(rep(y[1], hw), y, rep(y[n], hw))   # edges: replicate end points
  sapply(1:n, function(i) sum(cf * yp[i:(i + w - 1)]))
}
snv_hand <- function(y) (y - mean(y)) / sd(y)


## ===========================================================================
## 2. Replication: transcribed line by line from the CRAN sources; used as a
##    fallback when the packages are not installed
## ---------------------------------------------------------------------------
##  These are not "what I think the packages probably do". They were written
##  against
##    https://github.com/cran/baseline/blob/master/R/baseline.als.R
##    https://github.com/cran/prospectr/blob/master/R/savitzkyGolay.R
##  Every difference from the hand-written version is called out in a comment.
## ===========================================================================

## The core of baseline::baseline.als
##   Difference 1: lambda is a POWER OF TEN -- DD <- .create_DD(m) * 10^lambda
##   Difference 2: the weight update uses (y < z), not (y <= z), so a point
##                 sitting exactly on z gets weight 0
##   Difference 3: it breaks early once the weights stop changing (the
##                 hand-written version always runs the full niter)
##   Difference 4: the default is p = 0.05 (the hand-written version uses 0.01)
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

## The core of prospectr::savitzkyGolay
##   Difference 1: valid convolution -- the output is (w-1) points SHORTER than
##                 the input, (w-1)/2 removed from each end
##   Difference 2: it has a delta.wav argument; derivatives are divided by
##                 delta.wav^m so the units become "per cm-1"
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
## 3. Calling the real packages (only if installed)
## ===========================================================================

## ---- baseline ----------------------------------------------------------
## Three things to watch:
##   (a) spectra must be a matrix -- for a single spectrum use drop = FALSE,
##       otherwise baseline() stops with an error
##   (b) spectra go in ROWS, which matches the layout used throughout this repo
##   (c) lambda = 4 is what equals the hand-written lambda = 1e4
als_pkg <- function(Y, lambda = 4, p = 0.01, maxit = 20) {
  if (!is.matrix(Y)) Y <- matrix(Y, nrow = 1)
  bc <- baseline::baseline(Y, method = "als", lambda = lambda, p = p, maxit = maxit)
  baseline::getBaseline(bc)          # getCorrected(bc) returns the corrected spectra
}

## ---- prospectr ---------------------------------------------------------
## savitzkyGolay(X, m, p, w, delta.wav)
##   m = derivative order, p = polynomial order, w = window length (odd)
## standardNormalVariate(X) works row-wise with an n-1 denominator (same as sd())
sg_pkg  <- function(Y, m = 0, p = 3, w = 9, delta.wav = NULL) {
  if (is.null(delta.wav)) prospectr::savitzkyGolay(Y, m = m, p = p, w = w)
  else prospectr::savitzkyGolay(Y, m = m, p = p, w = w, delta.wav = delta.wav)
}
snv_pkg <- function(Y) prospectr::standardNormalVariate(Y)


## ===========================================================================
## 4. Item-by-item comparison
## ===========================================================================
line <- function(ch = "-") cat(strrep(ch, 74), "\n")
rmsd <- function(a, b) sqrt(mean((a - b)^2))
ok   <- function(v, tol = 1e-8) if (v < tol) "  <= identical" else ""

y  <- X[1, ]
gap <- 4                                    # (w-1)/2, w = 9

line("="); cat("A. ALS baseline\n"); line()

zh  <- als_hand(y)                          # lambda = 1e4
zr4 <- als_repl(y, lambda = 4, p = 0.01)    # replication, lambda = 4 -> 10^4
zr6 <- als_repl(y)                          # replication, package defaults
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
cat(sprintf("\n  1008 cm-1 peak height: hand %.1f | lambda=4 %.1f | lambda=6 (default) %.1f  -> %.1f%% apart\n",
            pk(zh), pk(zr4), pk(zr6), 100 * (pk(zh) - pk(zr6)) / pk(zh)))
cat("  >> lambda is log10. Pasting the hand-written 1e4 in here asks for 10^10000.\n")

line("="); cat("B. Savitzky-Golay\n"); line()

h  <- sg_hand(y)
r  <- sg_repl(y)
cat(sprintf("  length: hand %d | prospectr replication %d | input %d   (w-1 = 8 points lost)\n",
            length(h), length(r), length(y)))
cat(sprintf("  interior (drop %d each end) hand vs replication : max diff %.3e%s\n",
            gap, max(abs(h[(gap + 1):(length(y) - gap)] - r)),
            ok(max(abs(h[(gap + 1):(length(y) - gap)] - r)))))
if (HAS_PROSPECTR) {
  pp <- as.numeric(sg_pkg(y))
  cat(sprintf("  interior  hand vs prospectr::savitzkyGolay : max diff %.3e%s   <-- REAL PACKAGE\n",
              max(abs(h[(gap + 1):(length(y) - gap)] - pp)),
              ok(max(abs(h[(gap + 1):(length(y) - gap)] - pp)))))
}
cat(sprintf("  but the edges differ: hand point 1 = %.2f; prospectr has no such point.\n", h[1]))
cat("  >> after prospectr the wavenumber axis must be trimmed too:\n")
cat("     wn_sg <- wn[5:(length(wn)-4)]\n")

cat(sprintf("\n  1st-derivative scale: without delta.wav  max|d1| = %.4f  (per POINT)\n",
            max(abs(sg_repl(y, m = 1)))))
cat(sprintf("                        with delta.wav     max|d1| = %.4f  (per cm-1, a factor of %.4f)\n",
            max(abs(sg_repl(y, m = 1, delta.wav = STEP))), STEP))
cat("  >> Harmless for classification (a global constant), essential if you\n")
cat("     ever quote a derivative value.\n")

line("="); cat("C. SNV\n"); line()

A <- t(apply(X, 1, snv_hand))
B <- { C0 <- X - rowMeans(X); C0 / sqrt(rowSums(C0^2) / (ncol(X) - 1)) }
cat(sprintf("  hand vs prospectr replication : max diff %.3e%s\n", max(abs(A - B)), ok(max(abs(A - B)))))
if (HAS_PROSPECTR) {
  P <- snv_pkg(X)
  cat(sprintf("  hand vs standardNormalVariate : max diff %.3e%s   <-- REAL PACKAGE\n",
              max(abs(A - P)), ok(max(abs(A - P)))))
}
Cn <- (X - rowMeans(X)) / apply(X, 1, function(r) sqrt(mean((r - mean(r))^2)))
cat(sprintf("  had the denominator been n instead of n-1 : max diff %.3e  (factor sqrt(n/(n-1)) = %.6f)\n",
            max(abs(A - Cn)), sqrt(ncol(X) / (ncol(X) - 1))))
cat("  >> This one agrees exactly, because both use n-1.\n")


## ===========================================================================
## 5. What you would actually write: the whole chain, package version
## ===========================================================================
line("="); cat("D. The full preprocessing chain, package version\n"); line()

if (HAS_BASELINE && HAS_PROSPECTR) {
  bc    <- baseline::baseline(X, method = "als", lambda = 4, p = 0.01, maxit = 20)
  Xcorr <- baseline::getCorrected(bc)                       # 10 x 664
  Xsg   <- prospectr::savitzkyGolay(Xcorr, m = 0, p = 3, w = 9)   # 10 x 656
  wn_sg <- wn[5:(length(wn) - 4)]               # trim the axis too
  Xsnv  <- prospectr::standardNormalVariate(Xsg)
  Xd1   <- prospectr::savitzkyGolay(Xsnv, m = 1, p = 3, w = 9, delta.wav = STEP)
  wn_d1 <- wn_sg[5:(length(wn_sg) - 4)]               # trimmed AGAIN
  cat(sprintf("  Xcorr %d x %d -> Xsg %d x %d -> Xd1 %d x %d\n",
              nrow(Xcorr), ncol(Xcorr), nrow(Xsg), ncol(Xsg), nrow(Xd1), ncol(Xd1)))
  cat(sprintf("  wavenumber axis: %d -> %d -> %d\n", length(wn), length(wn_sg), length(wn_d1)))
  cat("  >> Every savitzkyGolay call costs 8 variables. A two-step chain costs 16.\n")
  cat("  >> When you map a VIP table back to wavenumbers, use the TRIMMED wn_d1,\n")
cat("     not the original wn.\n")
  write.csv(data.frame(wavenumber = wn_sg, t(Xsg)),
            "../data/rock_preprocessed_pkg.csv", row.names = FALSE)
  cat("  wrote ../data/rock_preprocessed_pkg.csv\n")
} else {
  cat("  This section runs only when both packages are installed.\n")
  cat("  install.packages(c(\"baseline\", \"prospectr\"))\n")
}

line("="); cat("Takeaways\n"); line()
cat("  1. lambda in the baseline package is log10 -- hand-written 1e4 is lambda = 4.\n")
cat("  2. prospectr::savitzkyGolay drops (w-1)/2 points at each end; trim the axis.\n")
cat("  3. SNV agrees exactly (both n-1) -- that one is a safe swap.\n")
cat("  4. Before mapping anything back to a chemical bond, check that your\n")
cat("     wavenumber vector and your data matrix are still aligned.\n")

## =====================================================================
## Data illustration: chemical process CCD (Montgomery, Design and Analysis
## of Experiments, Example 11.2). Coded factors: x1 = (time - 85)/5,
## x2 = (temp - 175)/5; axial distance 1.414, region = ball of radius sqrt(2).
## Output: ../results/data_example.rds and the complete LaTeX tables
##         ../tex/tables/chemdata.tex (Table 6), ../tex/tables/data.tex
##         (Table 7) and ../tex/tables/data_exact.tex (Table 8).
## Run from this folder: Rscript data_example.R
## =====================================================================
source("../R/slsecore.R")
## ---------------------------------------------------------- estimators
olse <- function(F, y) qr.coef(qr(F), y)

## Two-stage SLSE (Wang and Leblanc, 2008) with the optimal weight matrix
## W_i = V_i^{-1} evaluated at first-stage OLS estimates; Gauss-Newton fit.
slse <- function(F, y, maxit = 100) {
  n <- nrow(F); m <- ncol(F)
  th <- olse(F, y); r0 <- y - drop(F %*% th)
  s2 <- mean(r0^2); mu3 <- mean(r0^3); mu4 <- mean(r0^4)
  eta <- drop(F %*% th)
  v11 <- rep(s2, n); v12 <- mu3 + 2 * eta * s2
  v22 <- max(mu4 - s2^2, 1e-8) + 4 * eta^2 * s2 + 4 * eta * mu3
  det <- pmax(v11 * v22 - v12^2, 1e-10 * v11 * v22)
  w11 <- v22 / det; w12 <- -v12 / det; w22 <- v11 / det
  l11 <- sqrt(w11); l21 <- w12 / l11; l22 <- sqrt(pmax(w22 - l21^2, 1e-14))
  res <- function(g) {
    e <- drop(F %*% g[1:m]); r1 <- y - e; r2 <- y^2 - e^2 - g[m + 1]
    c(l11 * r1 + l21 * r2, l22 * r2)
  }
  jac <- function(g) {
    e <- drop(F %*% g[1:m])
    d1 <- -F; d2 <- -2 * e * F
    rbind(cbind(l11 * d1 + l21 * d2, -l21), cbind(l22 * d2, -l22))
  }
  g <- c(th, s2); Q <- sum(res(g)^2)
  for (it in 1:maxit) {
    J <- jac(g); step <- -qr.coef(qr(J), res(g)); step[is.na(step)] <- 0
    lam <- 1
    repeat {
      gn <- g + lam * step; Qn <- sum(res(gn)^2)
      if (Qn <= Q || lam < 1e-8) break
      lam <- lam / 2
    }
    conv <- abs(Q - Qn) <= 1e-12 * (1 + Q)
    g <- gn; Q <- Qn
    if (conv) break
  }
  g[1:m]
}

## ------------------------------------------------------- exact designs
efficient_round <- function(w, n) {          # Pukelsheim and Rieder (1992)
  l <- length(w); ni <- ceiling((n - l / 2) * w)
  while (sum(ni) < n) { j <- which.min(ni / w); ni[j] <- ni[j] + 1 }
  while (sum(ni) > n) { j <- which.max((ni - 1) / w); ni[j] <- ni[j] - 1 }
  ni
}

## ---------------------------------------------------------- data
d <- read.csv("../data/chemical_process.csv")
X <- as.matrix(d[, c("x1", "x2")]); y <- d$yield
F <- fvec(X)
b_ols <- olse(F, y); r <- y - drop(F %*% b_ols)
s2 <- mean(r^2); mu3 <- mean(r^3); mu4 <- mean(r^4)
that <- mu3^2 / (s2 * (mu4 - s2^2))
## parametric bootstrap from the fitted model with resampled residuals
set.seed(11)
tb <- replicate(2000, { e <- sample(r, length(r), TRUE); yy <- drop(F %*% b_ols) + e
  rr <- yy - drop(F %*% olse(F, yy)); a2 <- mean(rr^2); mean(rr^3)^2 / (a2 * (mean(rr^4) - a2^2)) })
cat("OLSE coef:", round(b_ols, 3), "\n")
cat(sprintf("sigma2=%.4f mu3=%.5f mu4=%.5f skew=%.3f t.hat=%.3f  boot 90%%: %.3f-%.3f\n",
            s2, mu3, mu4, mu3 / s2^1.5, that, quantile(tb, .05), quantile(tb, .95)))
b_sls <- slse(F, y); cat("SLSE coef:", round(b_sls, 3), "\n")

## the textbook CCD as an approximate design on the ball of radius sqrt(2):
## 5/13 at the centre (S0), 4/13 on the axial points (S1), 4/13 on the factorial points (Sk)
k <- 2; wccd <- c(5, 4, 4) / 13
Pm <- orbit_moments(k, "ball"); zccd <- drop(wccd %*% Pm)
f <- crit_fun("ball")
eff_ccd <- function(t) exp(rdesign_srs(k, t, "ball")$logpsi - f(zccd[1], zccd[2], k, t))
opt_hat <- rdesign_srs(k, that, "ball"); rob <- rdesign_srs(k, 0.725, "ball"); olsd <- rdesign_srs(k, 0, "ball")
effz <- function(z, t) exp(rdesign_srs(k, t, "ball")$logpsi - f(z[1], z[2], k, t))
tg <- seq(0, 0.95, by = 0.025)
out <- data.frame(design = c("Textbook CCD (5,4,4)", "R-opt SLSE, t = t.hat", "R-opt OLSE (t = 0)", "R-opt SLSE, t0 = 0.725"),
  w0 = c(wccd[1], opt_hat$weights[1], olsd$weights[1], rob$weights[1]),
  w1 = c(wccd[2], opt_hat$weights[2], olsd$weights[2], rob$weights[2]),
  wk = c(wccd[3], opt_hat$weights[3], olsd$weights[3], rob$weights[3]),
  eff_that = c(eff_ccd(that), 1, effz(olsd$moments, that), effz(rob$moments, that)),
  eff_05 = c(eff_ccd(.5), effz(opt_hat$moments, .5), effz(olsd$moments, .5), effz(rob$moments, .5)),
  eff_09 = c(eff_ccd(.9), effz(opt_hat$moments, .9), effz(olsd$moments, .9), effz(rob$moments, .9)),
  eff_min = c(min(sapply(tg, eff_ccd)), min(sapply(tg, function(t) effz(opt_hat$moments, t))),
              min(sapply(tg, function(t) effz(olsd$moments, t))), min(sapply(tg, function(t) effz(rob$moments, t)))))
print(out, digits = 4)
## 13-run exact versions (centre / axial / factorial counts)
for (dd in list(opt_hat, rob)) {
  ed <- expand_design(k, dd$orbits, dd$weights, "ball"); ni <- efficient_round(ed$w, 13)
  cat("13-run allocation (centre, each axial, each factorial):", ni, "\n")
  Xe <- ed$X[rep(1:nrow(ed$X), ni), ]
  cat("  eff at t.hat:", exp(rdesign_srs(k, that, "ball")$logpsi - logpsi_points(Xe, rep(1/13, 13), that)), "\n")
}
Xc <- rbind(matrix(0, 5, 2), orbit_points(2, 1, "ball"), orbit_points(2, 2, "ball"))
cat("textbook CCD exact eff at t.hat:", exp(rdesign_srs(k, that, "ball")$logpsi - logpsi_points(Xc, rep(1/13, 13), that)), "\n")
#saveRDS(list(that = that, tb = tb, out = out, b_ols = b_ols, b_sls = b_sls, mom = c(s2, mu3, mu4)), "../results/data_example.rds")

## symmetric exact designs: n0 centre runs, n1 replicates of each axial point,
## nk replicates of each factorial point, n0 + 4 n1 + 4 nk = N
sym_exact <- NULL
for (N in c(13, 17, 21, 25)) for (n1 in 1:6) for (nk in 1:6) {
  n0 <- N - 4 * n1 - 4 * nk; if (n0 < 1) next
  z <- drop(c(n0, 4 * n1, 4 * nk) / N) %*% orbit_moments(2, "ball")
  sym_exact <- rbind(sym_exact, data.frame(N, n0, n1, nk, eff_that = effz(z, that),
                                           eff_min = min(sapply(tg, function(t) effz(z, t)))))
}
best <- do.call(rbind, lapply(split(sym_exact, sym_exact$N), function(s) s[which.max(s$eff_that), ]))
bestr <- do.call(rbind, lapply(split(sym_exact, sym_exact$N), function(s) s[which.max(s$eff_min), ]))
print(sym_exact, digits = 4); print(best, digits = 4); print(bestr, digits = 4)
#res <- readRDS("../results/data_example.rds"); res$sym_exact <- sym_exact; saveRDS(res, "../results/data_example.rds")

## ---------------------------------------------------------- LaTeX tables
## Table 6: the data, in three column groups (factorial, centre, axial runs)
fx <- function(x) if (x == 0) "0" else sprintf("$%s$", format(x))
grp <- list(d[abs(d$x1) == 1 & abs(d$x2) == 1, ], d[d$x1 == 0 & d$x2 == 0, ],
            d[(d$x1 == 0) != (d$x2 == 0), ])
nr <- max(sapply(grp, nrow))
body <- sapply(seq_len(nr), function(i) paste(sapply(grp, function(g)
  if (i <= nrow(g)) sprintf("%s & %s & %.1f", fx(g$x1[i]), fx(g$x2[i]), g$yield[i]) else " & & "),
  collapse = " & "))

#dir.create("../tex/tables", showWarnings = FALSE)
writeLines(body, "../tables/chemdata.tex")
cat(body, sep = "\n")

## Table 7: orbit weights and R-efficiencies of the CCD and the R-optimal designs
labs <- c("Textbook CCD, weights $(5,4,4)/13$", "$\\xi^*_{\\SLS}(\\hat t)$",
          "$\\xi^*_{\\OLS}$", "$\\xi^*_{\\SLS}(0.725)$")
body <- sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\", labs, f4(out$w0), f4(out$w1), f4(out$wk),
                f4(out$eff_that), f4(out$eff_05), f4(out$eff_09), f4(out$eff_min))

#dir.create("../tex/tables", showWarnings = FALSE)
writeLines(body, "../tables/data.tex")
cat(body, sep = "\n")


## Table 8: for each run size N, the best symmetric exact design under each
## of two ranking criteria: (i) the R-efficiency at t = t.hat and (ii) the
## minimum R-efficiency over t in [0, 0.95] (maximin). If both criteria
## select the same allocation, it is shown once.
body <- c()
for (N in unique(sym_exact$N)) {
  r <- sym_exact[sym_exact$N == N, ]
  a <- r[which.max(r$eff_that), ]; b <- r[which.max(r$eff_min), ]
  same <- a$n0 == b$n0 && a$n1 == b$n1 && a$nk == b$nk
  rows <- if (same) list(list(a, "both")) else list(list(a, "R-eff at $\\hat t$"), list(b, "min over $t$"))
  for (i in seq_along(rows)) {
    x <- rows[[i]][[1]]
    body <- c(body, sprintf("%s & %s & (%d, %d, %d) & %s & %s \\\\", if (i == 1) N else "", rows[[i]][[2]],
                            x$n0, x$n1, x$nk, f4(x$eff_that), f4(x$eff_min)))
  }
  if (N != max(sym_exact$N)) body <- c(body, "\\addlinespace")
}

#dir.create("../tex/tables", showWarnings = FALSE)
writeLines(body, "../tables/data_exact.tex")
cat(body, sep = "\n")

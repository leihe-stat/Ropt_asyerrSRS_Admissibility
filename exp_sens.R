## =====================================================================
## Efficiency of designs under misspecified t and maximin choice of t0:
## Table 5 and Figures 2 and 3.
## Output: ../results/maximin.csv, ../results/sens.rds, the complete LaTeX
##         table ../tex/tables/maximin.tex, and ../figs/fig_olseff.pdf,
##         ../figs/fig_misspec.pdf. Run from this folder: Rscript exp_sens.R
## =====================================================================
source("../R/slsecore.R")
tg <- seq(0, 0.95, by = 0.025)
res <- NULL; Emats <- list()
for (region in c("cube", "ball")) for (k in 2:6) {
  f <- crit_fun(region)
  Z <- lapply(tg, function(t) rdesign_srs(k, t, region)$moments)
  opt <- mapply(function(z, t) f(z[1], z[2], k, t), Z, tg)
  E <- t(sapply(seq_along(tg), function(i) sapply(seq_along(tg), function(j) exp(opt[j] - f(Z[[i]][1], Z[[i]][2], k, tg[j])))))
  Emats[[paste(region, k)]] <- E
  minE <- apply(E, 1, min); i <- which.max(minE)
  res <- rbind(res, data.frame(region, k, ols_min = min(E[1, ]), t05_min = min(E[tg == 0.5, ]),
                               maximin_t0 = tg[i], maximin_eff = minE[i]))
}
print(res); write.csv(res, "../results/maximin.csv", row.names = FALSE)
saveRDS(list(tg = tg, E = Emats), "../results/sens.rds")

## Table 5 (LaTeX)
body <- c()
for (k in 2:6) {
  a <- res[res$region == "cube" & res$k == k, ]; b <- res[res$region == "ball" & res$k == k, ]
  body <- c(body, sprintf("%d & %s & %s & %.3f & %s & %s & %s & %.3f & %s \\\\", k,
                          f4(a$ols_min), f4(a$t05_min), a$maximin_t0, f4(a$maximin_eff),
                          f4(b$ols_min), f4(b$t05_min), b$maximin_t0, f4(b$maximin_eff)))
}
#dir.create("../tex/tables", showWarnings = FALSE)
writeLines(body, "../tables/maximin.tex")
cat(body, sep = "\n")

## Figure 2: R-efficiency of the OLSE-optimal design under the SLSE
tf <- seq(0, 0.99, by = 0.01)
pdf("../fig_olseff.pdf", width = 9, height = 4.2)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1))
for (region in c("cube", "ball")) {
  f <- crit_fun(region); cols <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
  plot(NA, xlim = c(0, 1), ylim = c(0.8, 1), xlab = expression(italic(t)), ylab = "R-efficiency",
       main = if (region == "cube") expression("(a) " * C[k]) else expression("(b) " * B[k]))
  grid()
  for (k in 2:6) {
    z0 <- rdesign_srs(k, 0, region)$moments
    ef <- sapply(tf, function(t) exp(rdesign_srs(k, t, region)$logpsi - f(z0[1], z0[2], k, t)))
    lines(tf, ef, col = cols[k - 1], lty = k - 1, lwd = 2)
  }
  legend("bottomleft", legend = paste("k =", 2:6), col = cols, lty = 1:5, lwd = 2, bty = "n")
}
dev.off()

## Figure 3: efficiency of xi*(t0) when the true value is t (k = 2 and 4, cube)
pdf("../fig_misspec.pdf", width = 9, height = 4.2)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1))
for (k in c(2, 4)) {
  E <- Emats[[paste("cube", k)]]
  plot(NA, xlim = c(0, 0.95), ylim = c(0.85, 1), xlab = expression("true " * italic(t)), ylab = "R-efficiency",
       main = bquote("(" * .(if (k == 2) "a" else "b") * ") " * C[.(k)]))
  grid()
  sel <- c(0, 0.25, 0.5, 0.75, 0.9); cols <- c("black", "#1f77b4", "#2ca02c", "#d62728", "#9467bd")
  for (s in seq_along(sel)) lines(tg, E[which.min(abs(tg - sel[s])), ], col = cols[s], lty = s, lwd = 2)
  legend("bottomleft", legend = as.expression(lapply(sel, function(v) bquote(italic(t)[0] ~ "=" ~ .(v)))),
         col = cols, lty = seq_along(sel), lwd = 2, bty = "n")
}
dev.off()

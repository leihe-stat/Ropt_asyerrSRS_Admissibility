## =====================================================================
## Table 4: CPU time of Algorithm 1 for t = 0.5 and k up to 40, compared
## with a CVX route to the same design: the reduced problem (18) is solved
## as a small semidefinite programme in CVXR (cvx_reduced.R) and the orbit
## weights with the fewest support points are then obtained by the same
## weight stage as in Algorithm 1 (min_support_weights).
## Both routes include the weight stage. CPU time = user + system time of
## the R process (the solvers run inside it): Algorithm 1 mean of 20 runs,
## CVX mean of 3 runs after one warm-up call.
## Output: ../results/compCVX.csv and the complete LaTeX table
##         ../tex/tables/compCVX.tex.
## Run from this folder: Rscript table5.R  (RERUN <- FALSE: table only)
## =====================================================================
source("/Users/wall/Desktop/2026-RSLSERespSur/0-EffCompRoptSRS/R/slsecore.R"); 
source("/Users/wall/Desktop/2026-RSLSERespSur/0-EffCompRoptSRS/R/cvx_reduced.R")
RERUN <- TRUE
t <- 0.5; ks <- c(2, 5, 10, 15, 20, 30, 40)
cpu <- function(expr, nrep) {
  t0 <- proc.time(); for (i in seq_len(nrep)) force(expr())
  sum((proc.time() - t0)[c("user.self", "sys.self")]) / nrep
}

if (RERUN) {
  invisible(rdesign_cvx(2, t, "cube", "SCS"))          # warm-up (loads CVXR and the solvers)
  for (i in 1:50) invisible(rdesign_srs(2, t))          # warm-up (byte compilation)
  sc <- NULL
  for (region in c("cube", "ball")) for (k in ks) {
    d  <- rdesign_srs(k, t, region)
    ta <- cpu(function() rdesign_srs(k, t, region), 20)
    res <- list()
    for (sv in c("SCS", "CLARABEL")) {
      dc <- rdesign_cvx(k, t, region, sv)
      res[[sv]] <- list(time = cpu(function() rdesign_cvx(k, t, region, sv), 3),
                        dlog = dc$logpsi - d$logpsi, status = dc$status,
                        same = identical(dc$orbits, d$orbits) && dc$npoints == d$npoints)
    }
    sc <- rbind(sc, data.frame(region, k, m = nparam(k), npoints = d$npoints,
                               orbits = paste(d$orbits, collapse = ","),
                               w1 = d$weights[1], w2 = d$weights[2], w3 = d$weights[3],
                               certificate = d$certificate, time = ta,
                               time_scs = res$SCS$time, dlog_scs = res$SCS$dlog,
                               status_scs = res$SCS$status, same_scs = res$SCS$same,
                               time_clarabel = res$CLARABEL$time, dlog_clarabel = res$CLARABEL$dlog,
                               status_clarabel = res$CLARABEL$status, same_clarabel = res$CLARABEL$same))
    print(tail(sc, 1)[, c("region", "k", "time", "time_scs", "time_clarabel", "dlog_scs", "dlog_clarabel", "same_scs", "same_clarabel")])
  }
  write.csv(sc, "/Users/wall/Desktop/2026-RSLSERespSur/0-EffCompRoptSRS/results/compCVX.csv", row.names = FALSE)
}

## ------------------------------------------------ LaTeX table
sc <- read.csv("/Users/wall/Desktop/2026-RSLSERespSur/0-EffCompRoptSRS/results/compCVX.csv")
np <- function(x) if (x > 1e9) { e <- floor(log10(x)); sprintf("$%.2f\\times10^{%d}$", x / 10^e, e) } else format(x, big.mark = ",")
ft <- function(x) if (x < 10) sprintf("%.3f", x) else sprintf("%.1f", x)
L <- c()
for (k in ks) {
  a <- sc[sc$region == "cube" & sc$k == k, ]; b <- sc[sc$region == "ball" & sc$k == k, ]
  L <- c(L, sprintf("%d & %d & %s & %s & %s & %s & %s & %s & %s & %s \\\\", k, a$m,
                    np(a$npoints), ft(a$time), ft(a$time_scs), ft(a$time_clarabel),
                    np(b$npoints), ft(b$time), ft(b$time_scs), ft(b$time_clarabel)))
}
#dir.create("../tex/tables", showWarnings = FALSE)
writeLines(L, "/Users/wall/Desktop/2026-RSLSERespSur/0-EffCompRoptSRS/tables/compCVX.tex")
cat(L, sep = "\n")

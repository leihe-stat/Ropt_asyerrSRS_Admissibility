## =====================================================================
## Robustness to starting values (Table 3).
## For every configuration (region, k, t), 100 starting values are drawn
## uniformly from the interior of the moment polygon (Z_C for the cube,
## Z_B for the ball). The SAME starting values are passed to
##   (i)  nleqslv, solving the equivalence-theorem equations directly, and
##   (ii) Algorithm 1 (rdesign_srs),
## so the two success rates are directly comparable, start by start.
## Output: ../results/robustness.csv and the complete LaTeX table
##         ../tex/tables/robustness.tex. Run from this folder: Rscript exp_algo.R
## =====================================================================
source("/Users/wall/Desktop/2026-RSLSERespSur/0-EffCompRoptSRS/R/slsecore.R"); library(nleqslv)
set.seed(2026)

n_start <- 100
tol     <- 1e-6          # tolerance on the moments, relative to the polygon scale

## uniform sample from the interior of the moment polygon (rejection sampling
## from its bounding box; the polygon is {u : U u > c} with unit normals)
draw_starts <- function(n, k, region) {
  P  <- orbit_moments(k, region)
  hc <- hull_constraints(P)
  lo <- apply(P, 2, min); hi <- apply(P, 2, max)
  out <- matrix(NA_real_, 0, 2)
  while (nrow(out) < n) {
    Z  <- cbind(runif(4 * n, lo[1], hi[1]), runif(4 * n, lo[2], hi[2]))
    ok <- apply(Z %*% t(hc$ui) - matrix(hc$ci, nrow(Z), length(hc$ci), byrow = TRUE) > 0, 1, all)
    out <- rbind(out, Z[ok, , drop = FALSE])
  }
  out[seq_len(n), , drop = FALSE]
}

## is z in D (positive-definite region), i.e. condition (12)?
in_D <- function(z, k, region) {
  if (region == "cube") { a <- z[1]; b <- z[2]
    b > 0 && b < a && a + (k - 1) * b - k * a^2 > 0
  } else { c <- z[1]; d <- z[2]; 0 < k * d && k * d < c && c < 1 / k }
}

## nleqslv from one start: returns "success", "no convergence",
## "root outside D" or "wrong root"
run_eq <- function(z0, k, t, region, ref, sc) {
  sysf <- eq_fun(region)
  r <- tryCatch(suppressWarnings(
         nleqslv(z0, function(u) sysf(u[1], u[2], k, t)$eq[2:3])),
       error = function(e) NULL)
  if (is.null(r) || r$termcd != 1 || any(!is.finite(r$x))) return("no convergence")
  if (all(abs((r$x - ref) / sc) < tol)) return("success")
  if (!in_D(r$x, k, region)) return("root outside D")
  "wrong root"
}

## Algorithm 1 from one start
run_alg <- function(z0, k, t, region, ref, sc) {
  z <- tryCatch(rdesign_srs(k, t, region, x0 = z0)$moments, error = function(e) c(NA, NA))
  if (all(is.finite(z)) && all(abs((z - ref) / sc) < tol)) "success" else "failure"
}

rob <- NULL
for (region in c("cube", "ball")) for (k in c(2, 3, 4, 6, 8, 11)) for (t in c(0, .5, .9)) {
  ref <- rdesign_srs(k, t, region)$moments
  sc  <- if (region == "cube") c(1, 1) else c(1 / k, 1 / k^2)
  Z0  <- draw_starts(n_start, k, region)
  eq  <- apply(Z0, 1, run_eq,  k = k, t = t, region = region, ref = ref, sc = sc)
  alg <- apply(Z0, 1, run_alg, k = k, t = t, region = region, ref = ref, sc = sc)
  rob <- rbind(rob, data.frame(
    region, k, t, n_start,
    eq_success       = mean(eq == "success"),
    proposed_success = mean(alg == "success"),
    eq_noconv        = mean(eq == "no convergence"),
    eq_outside_D     = mean(eq == "root outside D"),
    eq_wrong_root    = mean(eq == "wrong root"),
    both_success     = mean(eq == "success" & alg == "success")))
  print(tail(rob, 1))
}
write.csv(rob, "/Users/wall/Desktop/2026-RSLSERespSur/0-EffCompRoptSRS/results/startval.csv", row.names = FALSE)

## ------------------------------------------------ LaTeX table (Table 3)
## each cell: success rate (%) of nleqslv / of Algorithm 1
body <- c()
for (k in unique(rob$k)) {
  r <- rob[rob$k == k, ]
  cells <- c(sapply(c("cube", "ball"), function(g) sapply(c(0, .5, .9), function(t) {
    x <- r[r$region == g & r$t == t, ]
    sprintf("%.0f / %.0f", 100 * x$eq_success, 100 * x$proposed_success) })))
  body <- c(body, paste0(k, " & ", paste(cells, collapse = " & "), " \\\\"))
}

#dir.create("../tex/tables", showWarnings = FALSE)
writeLines(body, "/Users/wall/Desktop/2026-RSLSERespSur/0-EffCompRoptSRS/tables/startval.tex")
cat(body, sep = "\n")


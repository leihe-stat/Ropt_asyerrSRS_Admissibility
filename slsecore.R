## =====================================================================
## slsecore.R
## R-optimal designs under the second-order least squares estimator (SLSE)
## for second-order response surface (SRS) models on the cube C_k and the
## unit ball B_k, computed via design admissibility.
##
## The R-criterion is log psi = sum_j log e_j' M^{-1} e_j, M = G2 - t g1 g1'.
## Since M^{-1} = K B^{-1} K' and B is linear in the design, log psi is a
## convex function of the design moments (a, b) [cube] or (c, d) [ball].
##
## Algorithm 1: (1) Newton barrier method on the convex reduced problem;
## (2) minimal-support orbit weights; the residual of the equations of
## Theorem 3.3 is returned as an optimality certificate.
## =====================================================================

## ---------------------------------------------------------- regression
fvec <- function(X) {
  X <- as.matrix(X); if (ncol(X) == 1L && nrow(X) > 1L && FALSE) X <- t(X)
  k <- ncol(X); n <- nrow(X)
  out <- cbind(1, X^2, X)
  if (k > 1) for (i in 1:(k - 1)) for (j in (i + 1):k) out <- cbind(out, X[, i] * X[, j])
  unname(out)
}

nparam <- function(k) (k + 1) * (k + 2) / 2

info_mats <- function(F, w, t) {
  g  <- drop(crossprod(F, w))
  G2 <- crossprod(F * w, F)
  list(g = g, G2 = G2, M = G2 - t * tcrossprod(g))
}

logpsi_points <- function(X, w, t) {
  M <- info_mats(fvec(X), w, t)$M
  Mi <- tryCatch(solve(M), error = function(e) NULL)
  if (is.null(Mi)) return(Inf)
  d <- diag(Mi); if (any(d <= 0)) return(Inf)
  sum(log(d))
}

## sensitivity function phi(x, xi) of (A.1); xi is R-optimal iff max phi = m
sensitivity <- function(Xc, X, w, t) {
  im <- info_mats(fvec(X), w, t)
  Mi <- solve(im$M)
  A  <- Mi %*% diag(1 / diag(Mi)) %*% Mi
  Fc <- fvec(Xc)
  Fd <- sweep(Fc, 2, im$g)
  (1 - t) * rowSums((Fc %*% A) * Fc) + t * rowSums((Fd %*% A) * Fd)
}

## ------------------------------------------- closed-form reduced criteria
logpsi_cube <- function(a, b, k, t) {
  tau <- a + (k - 1) * b - k * a^2
  if (!(0 < b && b < a && a < 1 && tau > 0)) return(Inf)
  d0 <- (a - b + k * (b - a^2 * t)) / ((1 - t) * tau)
  dq <- (1 + (a^2 - b) / tau) / (a - b)
  if (d0 <= 0 || dq <= 0) return(Inf)
  log(d0) + k * log(dq) - k * log(a) - k * (k - 1) / 2 * log(b)
}

logpsi_ball <- function(c, d, k, t) {
  if (!(0 < k * d && k * d < c && c < 1 / k)) return(Inf)
  d0 <- (1 - k * c * t) / ((1 - t) * (1 - k * c))
  dq <- (1 + (c^2 - d) / (c * (1 - k * c))) / (c - k * d)
  if (dq <= 0) return(Inf)
  log(d0) + k * log(dq) - k * log(c) - k * (k - 1) / 2 * log(d)
}

crit_fun <- function(region) if (region == "cube") logpsi_cube else logpsi_ball

## ---------------------------------------------------------- orbits
orbit_moments <- function(k, region) {
  if (region == "cube") { l <- 0:k; cbind(l / k, l * (l - 1) / (k * (k - 1))) }
  else rbind(c(0, 0), c(1 / k, 0), c(1 / k, 1 / k^2))
}
orbit_size <- function(k, l, region) {     # l = orbit index starting at 0
  if (region == "cube") choose(k, l) * 2^l else c(1, 2 * k, 2^k)[l + 1]
}
orbit_points <- function(k, l, region) {
  signs <- function(r) as.matrix(expand.grid(rep(list(c(-1, 1)), r)))
  if (region == "cube") {
    if (l == 0) return(matrix(0, 1, k))
    S <- combn(k, l, simplify = FALSE); s <- signs(l)
    do.call(rbind, lapply(S, function(idx) { Z <- matrix(0, nrow(s), k); Z[, idx] <- s; Z }))
  } else {
    if (l == 0) return(matrix(0, 1, k))
    if (l == 1) return(rbind(diag(k), -diag(k)))
    unname(signs(k)) / sqrt(k)
  }
}

## ------------------------- Stage 1: Newton barrier method on the convex reduced problem
## The reduced criterion (13) [cube] or (14) [ball] is convex in the moments
## (Proposition 3.1). It is minimised over the moment polygon by a primal
## log-barrier (interior-point) method with damped Newton steps, using the
## exact gradient and Hessian obtained by symbolic differentiation (deriv()).
## Ball moments are rescaled as u = (k c, k^2 d), which maps the triangle
## onto the fixed triangle with vertices (0,0), (1,0), (1,1).
hull_constraints <- function(P) {
  h <- chull(P); h <- c(h, h[1]); ui <- NULL; ci <- NULL
  cen <- colMeans(P)
  for (i in seq_len(length(h) - 1)) {
    p <- P[h[i], ]; q <- P[h[i + 1], ]
    nrm <- c(q[2] - p[2], -(q[1] - p[1]))           # normal of edge p -> q
    if (sum(nrm * (cen - p)) < 0) nrm <- -nrm         # point inwards
    nrm <- nrm / sqrt(sum(nrm^2))
    ui <- rbind(ui, nrm); ci <- c(ci, sum(nrm * p))
  }
  list(ui = unname(ui), ci = ci)
}

## exact gradient and Hessian of the reduced criterion in scaled coordinates u
.crit_expr <- list(
  cube = quote(log((u1 - u2 + k * (u2 - u1^2 * t)) / ((1 - t) * (u1 + (k - 1) * u2 - k * u1^2))) +
               k * log((1 + (u1^2 - u2) / (u1 + (k - 1) * u2 - k * u1^2)) / (u1 - u2)) -
               k * log(u1) - k * (k - 1) / 2 * log(u2)),
  ball = quote(log((1 - u1 * t) / ((1 - t) * (1 - u1))) +
               k * log((1 + ((u1 / k)^2 - u2 / k^2) / ((u1 / k) * (1 - u1))) / (u1 / k - u2 / k)) -
               k * log(u1 / k) - k * (k - 1) / 2 * log(u2 / k^2)))
.crit_deriv <- lapply(.crit_expr, function(e) deriv(e, c("u1", "u2"), function(u1, u2, k, t) {}, hessian = TRUE))

crit_scaled <- function(u, k, t, region) {
  s <- if (region == "cube") c(1, 1) else c(1 / k, 1 / k^2)
  crit_fun(region)(u[1] * s[1], u[2] * s[2], k, t)
}
crit_ghess <- function(u, k, t, region) {
  v <- .crit_deriv[[region]](u[1], u[2], k, t)
  list(g = drop(attr(v, "gradient")), H = matrix(attr(v, "hessian"), 2, 2))
}

solve_reduced <- function(k, t, region = "cube", x0 = NULL,
                          mu0 = 1, shrink = 0.01, mu_min = 1e-12, tol = 1e-28, maxit = 50) {
  P  <- orbit_moments(k, region)
  s  <- if (region == "cube") c(1, 1) else c(1 / k, 1 / k^2)
  Pu <- sweep(P, 2, s, "/")
  hc <- hull_constraints(Pu); cen <- colMeans(Pu)
  f  <- function(u) crit_scaled(u, k, t, region)
  slack <- function(u) drop(hc$ui %*% u) - hc$ci
  u <- if (is.null(x0)) cen else x0 / s
  for (it in 1:60) {                     # make the start strictly feasible
    if (all(slack(u) > 1e-12) && is.finite(f(u))) break
    u <- 0.5 * u + 0.5 * cen
  }
  ## objective of the barrier subproblem: f(u) - mu * sum log slack(u)
  phi <- function(u, mu) { sl <- slack(u); if (any(sl <= 0)) return(Inf)
    v <- f(u); if (!is.finite(v)) return(Inf); v - mu * sum(log(sl)) }
  grad_phi <- function(u, mu) crit_ghess(u, k, t, region)$g - mu * colSums(hc$ui / slack(u))
  mu <- mu0; newton_steps <- 0
  repeat {
    for (it in 1:maxit) {                # damped Newton on phi(., mu)
      gh <- crit_ghess(u, k, t, region); sl <- slack(u)
      g <- gh$g - mu * colSums(hc$ui / sl)
      H <- gh$H + mu * crossprod(hc$ui / sl)
      ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
      if (min(ev) <= 1e-12 * max(abs(ev))) H <- H + (1e-12 * max(abs(ev)) - min(ev)) * diag(2)
      dlt <- -solve(H, g)
      lam2 <- -sum(g * dlt)              # squared Newton decrement
      if (lam2 / 2 <= tol) break
      ## backtracking: accept a step that satisfies the Armijo condition or,
      ## near the optimum where function values no longer resolve progress,
      ## one that reduces the gradient norm
      st <- 1; p0 <- phi(u, mu); gn0 <- sqrt(sum(g^2))
      while (st > 1e-16) {
        un <- u + st * dlt; pn <- phi(un, mu)
        if (is.finite(pn) && (pn <= p0 - 0.25 * st * lam2 ||
            sqrt(sum(grad_phi(un, mu)^2)) <= (1 - 0.25 * st) * gn0)) break
        st <- st / 2
      }
      if (st <= 1e-16) break
      u <- un; newton_steps <- newton_steps + 1
    }
    if (mu == 0) break
    mu <- if (mu * shrink < mu_min) 0 else mu * shrink   # finish with mu = 0
  }
  list(z = u * s, value = f(u), newton_steps = newton_steps,
       grad_norm = sqrt(sum(crit_ghess(u, k, t, region)$g^2)))
}

## ------------------------- equations (A.2) [cube] and (A.3) [ball]
h_system <- function(a, b, k, t) {
  tau <- a + (k - 1) * b - k * a^2
  lam <- a + (k - 2) * b - (k - 1) * a^2
  h1 <- -a / tau - a / lam - a * (a^2 - b) * k / (lam * tau)
  h2 <- tau / (lam * (a - b))
  h3 <- a^2 * (1 - t) / ((a - b + b * k - a^2 * k * t) * tau) +
        2 * (a^2 - b) / (lam * (a - b)) + (a^2 - b)^2 * k / (lam * tau * (a - b))
  H1 <- 2 * h1 + 1 / a - 2 * t * (h1 + a * k * h3 + a * h2)
  H0 <- -0.5 * a * k * (H1 + (k + 2) / a)
  list(eq = c(H0, h2 + h3 + H1, 2 * h3 + 1 / b), H1 = H1)
}

g_system <- function(c, d, k, t) {
  del <- c * (1 - k * c) + c^2 - d
  g1 <- -(del + c - k * d) / (del * (1 - k * c))
  g2 <- c * (1 - k * c) / (del * (c - k * d))
  g3 <- 2 * (c^2 - d) / (del * (c - k * d)) +
        k * (c^2 - d)^2 / (del * c * (1 - k * c) * (c - k * d)) +
        (1 - t) / ((1 - k * c) * (1 - k * c * t))
  G1 <- 2 * g1 + 1 / c - 2 * t * (g1 + c * g2 + k * c * g3)
  G0 <- -0.5 * c * k * (G1 + (k + 2) / c)
  list(eq = c(G0, g2 + g3 + G1, 2 * g3 + 2 * G1 + 1 / d), H1 = G1)
}

eq_fun <- function(region) if (region == "cube") h_system else g_system

## ------------- Stage 2: minimal-support orbit weights (Caratheodory / LP vertex)
min_support_weights <- function(k, z, region = "cube", tol = 1e-9) {
  ## Solve (17)/(18) for every subset of three orbits (a 3 x 3 Vandermonde-type
  ## system), keep the subsets with nonnegative weights, drop orbits whose
  ## weight is zero (this happens only if z lies on a chord p_i p_j), and
  ## return the design with the fewest support points.
  P <- orbit_moments(k, region); L <- nrow(P); rhs <- c(1, z); best <- NULL
  for (S in combn(L, 3, simplify = FALSE)) {
    A <- rbind(1, t(P[S, , drop = FALSE]))
    w <- tryCatch(solve(A, rhs), error = function(e) rep(NA, 3))
    if (anyNA(w) || any(w < -tol)) next
    keep <- w > tol
    S <- S[keep]; w <- w[keep]
    np <- sum(sapply(S - 1, function(l) orbit_size(k, l, region)))
    better <- is.null(best) || np < best$npoints ||
      (np == best$npoints && sum(S) > sum(best$orbits + 1))
    if (better) best <- list(npoints = np, orbits = S - 1, weights = w / sum(w))
  }
  best
}

expand_design <- function(k, orbits, weights, region = "cube") {
  X <- NULL; w <- NULL
  for (i in seq_along(orbits)) {
    Q <- orbit_points(k, orbits[i], region)
    X <- rbind(X, Q); w <- c(w, rep(weights[i] / nrow(Q), nrow(Q)))
  }
  list(X = X, w = w)
}

## ------------------------------------------------- Algorithm 1 (user level)
## Stage 1: Newton barrier method for the moments; Stage 2: minimal-support
## weights; the residuals of the equations of Theorem 3.3 at the returned
## moments are reported as an optimality certificate.
rdesign_srs <- function(k, t, region = c("cube", "ball"), x0 = NULL) {
  region <- match.arg(region)
  st1 <- solve_reduced(k, t, region, x0)
  z <- st1$z
  ms <- min_support_weights(k, z, region)
  cert <- max(abs(eq_fun(region)(z[1], z[2], k, t)$eq))
  structure(list(k = k, t = t, region = region, moments = z, logpsi = st1$value,
                 orbits = ms$orbits, weights = ms$weights, npoints = ms$npoints,
                 newton_steps = st1$newton_steps, grad_norm = st1$grad_norm,
                 certificate = cert),
            class = "srsdesign")
}

print.srsdesign <- function(x, ...) {
  lab <- if (x$region == "cube") paste0("E", x$orbits) else c("S0", "S1", "Sk")[x$orbits + 1]
  cat(sprintf("R-optimal SLSE design, %s, k = %d, t = %.3f (%d support points)\n",
              x$region, x$k, x$t, x$npoints))
  print(setNames(round(x$weights, 4), lab))
  invisible(x)
}

## R-efficiency (paper's definition: ratio of criterion values)
reff <- function(design, t) {
  f <- crit_fun(design$region)
  opt <- rdesign_srs(design$k, t, design$region)
  exp(opt$logpsi - f(design$moments[1], design$moments[2], design$k, t))
}

## ------------------------------------------------------ LaTeX tables
## Write a complete LaTeX table (table environment with caption, label,
## booktabs header and body) to `file`, ready for \input in the manuscript.
## header and body are character vectors of table rows (without \toprule
## etc.); `pre` holds optional lines placed before the tabular (e.g. a
## \setlength). Requires the LaTeX packages booktabs and adjustbox.
latex_table <- function(file, body, header, colspec, caption, label,
                        pos = "t", size = "\\footnotesize", pre = NULL) {
  L <- c(sprintf("\\begin{table}[%s]", pos), "\\centering",
         sprintf("\\caption{%s}\\label{%s}", caption, label), size, pre,
         "\\begin{adjustbox}{max width=\\textwidth}",
         sprintf("\\begin{tabular}{%s}", colspec), "\\toprule", header, "\\midrule",
         body, "\\bottomrule", "\\end{tabular}", "\\end{adjustbox}", "\\end{table}")
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  writeLines(L, file)
  invisible(L)
}
f4 <- function(x) formatC(x, format = "f", digits = 4)

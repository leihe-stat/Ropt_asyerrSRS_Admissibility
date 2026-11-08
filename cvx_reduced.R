## =====================================================================
## The reduced problem (18) solved by CVX (CVXR) as the SDP (A.6) of
## Proposition A.3, followed by the same
## minimal-support weight stage as Algorithm 1. Used for Table 4.
##
## For an invariant design with moments z, the matrix B of (7) is block
## diagonal (after reordering) with blocks
##   * the (k+2) x (k+2) block B_q(z) for the augmenting coordinate, the
##     intercept and the squares x_1^2, ..., x_k^2,
##   * k scalar blocks E x_i^2 for the linear terms,
##   * k(k-1)/2 scalar blocks E x_i^2 x_j^2 for the interactions.
## Hence, with e_j the j-th unit vector,
##   log psi = log(e_2' B_q^{-1} e_2) + k log(e_3' B_q^{-1} e_3)
##             - k log(E x_i^2) - k(k-1)/2 log(E x_i^2 x_j^2),
## and, by the Schur complement (Lemma A.2), minimising log psi over the
## polygon is the convex programme
##   max  log r0 + k log r2 + k log m2 + k(k-1)/2 log m22
##   s.t. B_q(z) - r0 e_2 e_2' >= 0,  B_q(z) - r2 e_3 e_3' >= 0,  U u >= c,
## where B_q is affine in the moments. The entries of B_q are
##   B_q = [1, sqrt(t) g'; sqrt(t) g, G],  g = (1, m2 1_k'),
##   G   = [1, m2 1_k'; m2 1_k, (m4 - m22) I_k + m22 1_k 1_k'],
## with m2 = E x_i^2, m22 = E x_i^2 x_j^2 and m4 = E x_i^4; on the cube
## (orbits E_l) m4 = m2 = a, m22 = b; on the ball (orbits S_0, S_1, S_k)
## m2 = c, m22 = d and m4 = c - (k-1) d.
## =====================================================================
suppressMessages(library(CVXR))

cvx_reduced <- function(k, t, region = c("cube", "ball"), solver = "CLARABEL") {
  region <- match.arg(region)
  s  <- if (region == "cube") c(1, 1) else c(1 / k, 1 / k^2)   # u = z / s (scaled moments)
  hc <- hull_constraints(sweep(orbit_moments(k, region), 2, s, "/"))
  n  <- k + 2; sq <- 3:n
  A0 <- matrix(0, n, n); A0[1, 1] <- 1; A0[1, 2] <- A0[2, 1] <- sqrt(t); A0[2, 2] <- 1
  A2 <- matrix(0, n, n)                                         # coefficient of m2
  A2[1, sq] <- A2[sq, 1] <- sqrt(t); A2[2, sq] <- A2[sq, 2] <- 1
  A22 <- matrix(0, n, n); A22[sq, sq] <- 1                      # coefficient of m22
  D <- cbind(sq, sq)                                            # diagonal of the squares block
  A2[D]  <- 1                                                   # m4 = m2 + (cube 0 / ball -(k-1)) m22
  A22[D] <- if (region == "cube") 0 else -(k - 1)
  E2 <- matrix(0, n, n); E2[2, 2] <- 1
  E3 <- matrix(0, n, n); E3[3, 3] <- 1

  u  <- Variable(2); r0 <- Variable(); r2 <- Variable()
  m2 <- u[1] * s[1]; m22 <- u[2] * s[2]
  Bq <- A0 + m2 * A2 + m22 * A22
  obj <- log(r0) + k * log(r2) + k * log(u[1]) + k * (k - 1) / 2 * log(u[2])
  cons <- list(hc$ui %*% u >= hc$ci, (Bq - r0 * E2) %>>% 0, (Bq - r2 * E3) %>>% 0)
  prob <- Problem(Maximize(obj), cons)
  res  <- psolve(prob, solver = solver)
  z <- drop(value(u)) * s
  list(z = z, status = status(prob), logpsi = crit_fun(region)(z[1], z[2], k, t))
}

## CVX route to the design: (18) by CVX, then the weight stage of Algorithm 1
rdesign_cvx <- function(k, t, region = c("cube", "ball"), solver = "CLARABEL") {
  region <- match.arg(region)
  st <- cvx_reduced(k, t, region, solver)
  ms <- min_support_weights(k, st$z, region)
  list(k = k, t = t, region = region, moments = st$z, logpsi = st$logpsi, status = st$status,
       orbits = ms$orbits, weights = ms$weights, npoints = ms$npoints)
}

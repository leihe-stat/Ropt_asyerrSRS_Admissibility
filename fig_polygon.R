## =====================================================================
## Figure 1: the moment polygon Z_C and the positive-definiteness region
## (12) for the SRS model on the cube, k = 3.
##   (a) full view;  (b) enlargement of the lower boundary.
## Self-contained: needs base R only. Output: fig_polygon.pdf
## =====================================================================

k  <- 3
tt <- 0.5                                   # skewness parameter for the marked design

## vertices p_l = (l/k, l(l-1)/(k(k-1))), l = 0, ..., k
l <- 0:k
P <- cbind(a = l / k, b = l * (l - 1) / (k * (k - 1)))

## parabola tau = a + (k-1) b - k a^2 = 0, i.e. b = (k a^2 - a)/(k - 1)
par_b <- function(a) (k * a^2 - a) / (k - 1)

## closed-form reduced R-criterion (13) on the cube
logpsi_cube <- function(a, b, k, t) {
  tau <- a + (k - 1) * b - k * a^2
  if (!(0 < b && b < a && a < 1 && tau > 0)) return(Inf)
  d0 <- (a - b + k * (b - a^2 * t)) / ((1 - t) * tau)
  dq <- (1 + (a^2 - b) / tau) / (a - b)
  log(d0) + k * log(dq) - k * log(a) - k * (k - 1) / 2 * log(b)
}

## R-optimal moments for t = tt: minimise (13) over Z_C, parametrising the
## orbit weights by a softmax so that the moments stay in the polygon
obj <- function(v) {
  w <- exp(c(0, v)); w <- w / sum(w)
  z <- drop(w %*% P)
  logpsi_cube(z[1], z[2], k, tt)
}
fit <- optim(rep(0, k), obj, method = "BFGS", control = list(reltol = 1e-14, maxit = 5000))
w <- exp(c(0, fit$par)); w <- w / sum(w)
z <- drop(w %*% P)                          # (a*, b*) = (0.7362, 0.5828)

## colours (fills light, borders saturated) and drawing helpers
blue_f <- "#cde2fb"; blue <- "#2a78d6"
orng_f <- "#f9d3c3"; orng <- "#eb6834"
ink <- "#0b0b0b"; grey <- "#52514e"

panel <- function(xlim, ylim, main, full) {
  plot(NA, xlim = xlim, ylim = ylim, xlab = "a", ylab = "b", las = 1, main = main,
       font.main = 1, cex.main = 1, asp = if (full) 1 else NA)
  ## strips that satisfy (12) but lie outside Z_C: between tau = 0 and the
  ## lower edges p_l p_{l+1}, l = 1, ..., k - 1 (hatched for print/CVD)
  for (j in 1:(k - 1)) {
    aa <- seq(P[j + 1, 1], P[j + 2, 1], length.out = 300)
    chord <- approx(P[(j + 1):(j + 2), 1], P[(j + 1):(j + 2), 2], aa)$y
    xs <- c(aa, rev(aa)); ys <- c(pmax(par_b(aa), 0), rev(chord))
    polygon(xs, ys, col = orng_f, border = NA)
    polygon(xs, ys, density = if (full) 40 else 20, angle = 45, col = orng,
            border = NA, lwd = 0.8)
  }
  polygon(P[, 1], P[, 2], col = blue_f, border = blue, lwd = 2)     # Z_C
  aa <- seq(1 / k, 1, length.out = 400)
  lines(aa, par_b(aa), lwd = 2, col = ink)                          # tau = 0
  abline(0, 1, lty = 2, col = grey)                                 # b = a
  points(P[, 1], P[, 2], pch = 19, cex = 1.1)
  lab <- sapply(l, function(i) as.expression(bquote(bold(p)[.(i)])))
  text(P[, 1], P[, 2], labels = lab, pos = c(1, 1, 4, 2), offset = 0.6, xpd = TRUE)
}

pdf("../fig_polygon.pdf", width = 7.4, height = 4.2)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1), oma = c(2.2, 0, 0, 0))

## (a) full view
panel(c(0, 1), c(-0.06, 1), "(a)", TRUE)
points(z[1], z[2], pch = 23, bg = "white", cex = 1.4, lwd = 1.5)
segments(0.50, 0.70, z[1] - 0.015, z[2] + 0.01, col = grey)
text(0.50, 0.70, expression(xi[SLS]^"*" * (0.5)), pos = 2, offset = 0.2)
text(0.95, 0.68, expression(tau == 0), cex = 0.95)
text(0.16, 0.24, "b = a", srt = 45, col = grey, cex = 0.9)
rect(0.30, -0.02, 0.72, 0.40, lty = 3, border = grey)

## (b) enlargement of the dotted box
panel(c(0.30, 0.72), c(-0.02, 0.40), "(b) enlargement of the dotted box in (a)", FALSE)
text(0.58, 0.155, expression(tau == 0), cex = 0.95, pos = 4)
text(0.47, 0.20, expression("chord " * bold(p)[1] * bold(p)[2]), cex = 0.9, col = blue)

## shared legend below both panels
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", horiz = TRUE, bty = "n", cex = 0.9, xpd = TRUE,
       legend = c(expression(Z[C] * " (moments of invariant designs)"),
                  expression("satisfies (13) but lies outside " * Z[C])),
       fill = c(blue_f, orng), density = c(NA, 40), angle = c(0, 45),
       border = c(blue, orng))
dev.off()

cat(sprintf("R-optimal moments for t = %.1f: a = %.6f, b = %.6f\n", tt, z[1], z[2]))

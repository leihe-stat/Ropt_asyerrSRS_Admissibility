## =====================================================================
## Tables 1 and 2: R-optimal SLSE designs with minimal support on the cube
## (Table 1) and the ball (Table 2) for 2 <= k <= 11 and t = 0, 0.5, 0.9.
## Output: ../results/tables.csv and the complete LaTeX tables
##         ../tex/tables/designs_cube.tex and ../tex/tables/designs_ball.tex.
## Run from this folder: Rscript exp_scal.R
## =====================================================================
source("../R/slsecore.R")
tab <- NULL
for (region in c("cube", "ball")) for (k in 2:11) for (t in c(0, .5, .9)) {
  d <- rdesign_srs(k, t, region)
  tab <- rbind(tab, data.frame(region, k, t, orbits = paste(d$orbits, collapse = ","),
                               w1 = d$weights[1], w2 = d$weights[2], w3 = d$weights[3], npoints = d$npoints,
                               m1 = d$moments[1], m2 = d$moments[2]))
}
write.csv(tab, "../results/tables1_2.csv", row.names = FALSE)

## ------------------------------------------------ LaTeX tables
## two blocks side by side: k = 2..6 on the left, k = 7..11 on the right
design_rows <- function(region) {
  s <- tab[tab$region == region, ]; ks <- 2:11; L <- c()
  half <- function(k) {
    r <- s[s$k == k, ]
    sapply(1:3, function(i) {
      lead <- if (i == 2) as.character(k) else ""
      if (region == "cube")
        sprintf("%s & %s & (%s) & %s & %s & %s", lead, format(r$t[i]), gsub(",", ", ", r$orbits[i]),
                f4(r$w1[i]), f4(r$w2[i]), f4(r$w3[i]))
      else sprintf("%s & %s & %s & %s & %s", lead, format(r$t[i]), f4(r$w1[i]), f4(r$w2[i]), f4(r$w3[i]))
    })
  }
  for (j in 1:5) {
    L <- c(L, paste0(half(ks[j]), " & ", half(ks[j + 5]), " \\\\"))
    if (j < 5) L <- c(L, "\\midrule")
  }
  L
}
#dir.create("../tex/tables", showWarnings = FALSE)
writeLines(design_rows("cube"), "../tables/designs_cube.tex")
cat(design_rows("cube"), sep = "\n")

writeLines(design_rows("ball"), "../tables/designs_ball.tex")
cat(design_rows("ball"), sep = "\n")

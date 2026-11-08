# R-opt_asyerrSRS_Admissibility

The code computes R-optimal approximate designs under the second-order least squares
estimator (SLSE) for second-order response surface (SRS) models in `k` factors on the
cube `[-1, 1]^k` and the unit ball. By design admissibility, the problem reduces to the
minimization of a convex function of two moments over a polygon. **Algorithm 1** of the
paper solves this reduced problem by a Newton interior-point (barrier) method and then
recovers the design with the fewest support points from a small linear system.

## Requirements

| Software | Needed for |
|---|---|
| R (≥ 4.0; the paper used R 4.6.1) | all scripts |
| [`nleqslv`](https://CRAN.R-project.org/package=nleqslv) | `exp_startval.R` (Table 3) |
| [`CVXR`](https://CRAN.R-project.org/package=CVXR) (≥ 1.9) with the solvers SCS and Clarabel | `cvx_reduced.R`, `exp_compareCVX.R` (Table 4) |

Everything else uses base R only.

```r
install.packages(c("nleqslv", "CVXR"))
```

## Repository structure

```
R/
  slsecore.R          core functions: reduced criteria, Algorithm 1, sensitivity function
  cvx_reduced.R       CVX method: semidefinite programme (24) solved with CVXR
  fig_polygon.R       Figure 1   (moment polygon for k = 3; base R only)
  exp_algo.R          Tables 1-2 (R-optimal designs on the cube and the ball, 2 <= k <= 11)
  exp_startval.R      Table 3    (robustness to starting values: Algorithm 1 vs nleqslv)
  exp_compareCVX.R    Table 4    (CPU time: Algorithm 1 vs the CVX method, k up to 40)
  exp_sens.R          Table 5, Figures 2-3 (OLSE design efficiency, misspecified t, maximin t0)
  data_example.R      Tables 6-8 (chemical process example, Section 6)
data/
  chemical_process.csv   Montgomery (2017), Example 11.2: columns x1, x2, yield (13 runs)
results/              CSV/RDS output of the scripts
tables/               LaTeX table bodies written by the scripts
```

## Quick start

```r
source("R/slsecore.R")

# R-optimal SLSE design on the cube, k = 3 factors, t = 0.5
d <- rdesign_srs(k = 3, t = 0.5, region = "cube")
d
#> R-optimal SLSE design, cube, k = 3, t = 0.500 (21 support points)
#>     E0     E2     E3
#> 0.1104 0.4603 0.4294

d$moments       # optimal moments (a, b) = (0.736224, 0.582805)
d$certificate   # residual of the equations of Theorem A.1 (about 1e-15)

# on the unit ball, k = 2
rdesign_srs(k = 2, t = 0.5, region = "ball")
#> R-optimal SLSE design, ball, k = 2, t = 0.500 (9 support points)
#>     S0     S1     Sk
#> 0.2849 0.3147 0.4004

# R-efficiency of the OLSE-optimal design (t = 0) when the true t is 0.5
reff(rdesign_srs(3, 0, "cube"), t = 0.5)
#> 0.9963

# support points and weights of the design
expand_design(3, d$orbits, d$weights, "cube")
```

Orbits are labelled as in Section 3 of the paper. On the cube, `E_l` is the set of points
with `l` coordinates equal to ±1 and the others equal to 0. On the ball, `S0` is the
centre, `S1` consists of the points ±e_i, and `Sk` consists of the points with all
coordinates equal to ±1/√k.

## Main functions (`slsecore.R`)

| Function | Description |
|---|---|
| `rdesign_srs(k, t, region, x0 = NULL)` | Algorithm 1: optimal moments, orbit weights, number of support points, Newton steps, certificate |
| `solve_reduced(k, t, region, x0)` | Newton interior-point stage: barrier method on the reduced problem (17) |
| `min_support_weights(k, z, region)` | Weight stage: solves (18) or (19) for all triples of orbits and returns the design with the fewest support points |
| `logpsi_cube(a, b, k, t)`, `logpsi_ball(c, d, k, t)` | Reduced R-criteria (14) and (15) |
| `h_system`, `g_system` | Equations (A.2) and (A.3) of Theorem A.1, used for verification |
| `sensitivity(Xc, X, w, t)` | Sensitivity function (A.1) of the equivalence theorem |
| `logpsi_points(X, w, t)` | log R-criterion of an arbitrary design given by points `X` and weights `w` |
| `reff(design, t)` | R-efficiency of a design under the SLSE with parameter `t` |
| `orbit_points`, `orbit_moments`, `orbit_size`, `expand_design` | Orbits and their moments, and expansion of an invariant design into its support points |

`cvx_reduced.R` provides `rdesign_cvx(k, t, region, solver)`, which solves the reduced
problem as the semidefinite programme (24) in CVXR (`solver = "SCS"` or `"CLARABEL"`) and
then applies the same weight stage.

## Reproducing the results

The scripts source `slsecore.R` (and `cvx_reduced.R`) and read or write the folders
`data/`, `results/` and `tables/`. **Before running them, set the paths at the top of each
script to the location of your copy of the repository.** For example, replace the absolute
paths with paths relative to the `R/` folder:

```r
source("slsecore.R")
d <- read.csv("../data/chemical_process.csv")
write.csv(tab, "../results/tables1_2.csv", row.names = FALSE)
```

Then run the scripts from the `R/` folder:

| Script | Reproduces | Approximate run time |
|---|---|---|
| `fig_polygon.R` | Figure 1 | seconds |
| `exp_algo.R` | Tables 1–2 | seconds |
| `exp_startval.R` | Table 3 | < 1 min |
| `exp_compareCVX.R` | Table 4 | about 10 min |
| `exp_sens.R` | Table 5, Figures 2–3 | a few minutes |
| `data_example.R` | Tables 6–8 | < 1 min |

The scripts write the body of each LaTeX table to `tables/` (for `\input` inside a
`tabular` environment) and the underlying numbers to `results/`. Random starting values
(Table 3) and the bootstrap interval for `t` (Section 6) use fixed seeds. CPU times
(Table 4) depend on the machine. The times in the paper were obtained with R 4.6.1 on
macOS with 8 GB of RAM.

## Data

`data/chemical_process.csv` contains the chemical process experiment of Montgomery (2017,
*Design and Analysis of Experiments*, 9th ed., Example 11.2). It is a rotatable central
composite design with 13 runs, given in the coded factors
`x1 = (time - 85)/5` and `x2 = (temperature - 175)/5`, with the process yield as response.

## Contact

Lei He, Department of Statistics, Anhui Normal University, Wuhu, China
(lhstat@163.com).

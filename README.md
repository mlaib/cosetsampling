# cosetsampling

[![License: GPL (>= 3)](https://img.shields.io/badge/license-GPL%20(%3E%3D%203)-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R >= 4.0.0](https://img.shields.io/badge/R-%3E%3D%204.0.0-276DC3.svg)](https://www.r-project.org/)
[![Version 0.1.0](https://img.shields.io/badge/version-0.1.0-informational.svg)](DESCRIPTION)

**Exactly balanced sampling on overlapping stratifications with coset designs.**

`cosetsampling` implements *coset designs*. These are probability sampling designs that are exactly balanced on
several overlapping stratifications at once, with first- and second-order inclusion probabilities in closed form.
They apply when the units of the frame can be indexed by a vector space over a prime field and the
stratifications are coset partitions. Examples are complete factorial frames stratified by their margins,
cross-classifications or linear combinations of the factors, and square grids of prime side stratified by rows,
columns and wrapped diagonals.

A coset design draws an admissible direction, a subspace transversal to every stratification, and then `r` random
cosets of it. Every sample therefore meets each stratum in exactly its proportional share. The package also
evaluates samples from other designs, such as the cube method, on the same scale.

## Installation

Install the development version from GitHub with [remotes](https://cran.r-project.org/package=remotes):

```r
# install.packages("remotes")
remotes::install_github("mlaib/cosetsampling")
```

To build the vignette as well, which requires knitr and rmarkdown:

```r
remotes::install_github("mlaib/cosetsampling", build_vignettes = TRUE)
```

The package depends only on base R (>= 4.0.0). `minimax_mixture()` uses the optional package
[lpSolve](https://cran.r-project.org/package=lpSolve): `install.packages("lpSolve")`.

## Quick start

A 61 × 61 grid (61 is prime), stratified by rows, columns and both families of wrapped diagonals:

```r
library(cosetsampling)

g <- grid_frame(61)
S <- strata_grid(g, c("rows", "columns", "antidiagonals", "diagonals"))

# Equal-weight coset design: one random admissible line, one random coset (k = 61 units)
d <- coset_design(g, S, h = 1)

# Exact design effect for a response, here the central block of R's volcano data
y <- as.vector(t(volcano[14:74, ]))
design_effect(d, y)
#> [1] 0.08782422

# Draw a sample; it is exactly balanced on all four stratifications
s <- draw_sample(d)
balance_check(g, S, s$units)$exact
#> [1] TRUE
```

With `r >= 2` cosets the design admits an unbiased variance estimator:

```r
d3 <- coset_design(g, S, h = 1, r = 3)   # k = 183
s3 <- draw_sample(d3)
random_group_variance(s3, y, d3)         # estimate, variance, degrees of freedom and interval
```

Complete factorial frames work the same way, and `minimax_mixture()` finds the mixture of directions with the
smallest worst-case variance:

```r
f <- coset_frame(2, 7)                   # the 2^7 factorial frame
M <- strata_margins(f, order = 1)        # balance on all seven one-variable margins
m <- minimax_mixture(f, M, h = 3)        # samples of size 2^3 = 8; needs lpSolve
m$ratio                                  # worst case relative to the lower bound
#> [1] 1
```

## Main functions

| Task | Functions |
|---|---|
| Frames and stratifications | `grid_frame()`, `coset_frame()`, `factorial_frame()`, `strata_grid()`, `strata_margins()`, `strata_linear()` |
| Admissible directions | `admissible_exists()`, `admissible_directions()`, `sample_direction()`, `is_admissible()` |
| Designs and samples | `coset_design()`, `minimax_mixture()`, `draw_sample()` |
| Inclusion probabilities and variances | `inclusion_probabilities()`, `joint_inclusion()`, `spectral_weights()`, `design_variance()`, `design_effect()`, `srs_variance()` |
| Robustness | `worst_case_factor()`, `is_equal_weight()` |
| Variance estimation | `random_group_variance()`, `regression_variance()` |
| Samples from any design, for example the cube method | `balance_check()`, `main_effect_weight()`, `main_effect_share()`, `worst_case_bound()` |
| Finite-field utilities | `rref_mod()`, `rank_mod()`, `nullspace_mod()`, `gaussian_binomial()`, `stratum_labels()`, `main_effect_set()`, `main_effect_count()`, `fourier_power()` |

## Documentation

- Reference manual: `help(package = "cosetsampling")`.
- Vignette: `vignette("cosetsampling")`, if built at installation.
- Worked examples of every function: [`notebooks/cosetsampling_functions.ipynb`](notebooks/cosetsampling_functions.ipynb),
  also as a stand-alone page, [`notebooks/cosetsampling_functions.html`](notebooks/cosetsampling_functions.html).

## Scope

- The field order `q` must be prime. A frame over a prime-power field can be written over its prime field, with
  its stratifications given as linear forms through `strata_linear()`.
- Sample sizes are `r * q^h`, and every unit has inclusion probability `k / N`.
- Exact balance holds for stratifications defined by cosets of subspaces; other auxiliary information is outside
  the scope of the package.

## Reference

The theory behind the package is developed in:

> Laib, M. (2026). *Exactly balanced sampling on overlapping stratifications with affine flats*. Manuscript.

The code that reproduces every result of the manuscript is in
[mlaib/balanced-sampling-affine-flats](https://github.com/mlaib/balanced-sampling-affine-flats).

## Licence

GPL (>= 3). See [LICENSE.md](LICENSE.md).

## Author

Mohamed Laib, Luxembourg Institute of Science and Technology (LIST).
Contact: laib.med@gmail.com, mohamed.laib@list.lu. Bug reports and questions are welcome as
[GitHub issues](https://github.com/mlaib/cosetsampling/issues).

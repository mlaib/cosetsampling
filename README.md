# cosetsampling

Exactly balanced sampling on overlapping stratifications with coset designs.

When the units of a frame can be indexed by a vector space over a prime field, as in complete factorial frames
and square grids of prime side, and the stratifications are coset partitions (margins, cross-classifications,
rows, columns, wrapped diagonals), a coset design selects samples that are exactly balanced on all the
stratifications at once, with closed-form first- and second-order inclusion probabilities. The package:

- lists admissible directions, or draws them by rejection (`admissible_directions()`, `sample_direction()`);
- builds coset designs with uniform, user-given or minimax mixtures of directions (`coset_design()`,
  `minimax_mixture()`) and draws samples (`draw_sample()`);
- gives inclusion probabilities, spectral weights, exact variances, design effects and worst-case factors
  (`joint_inclusion()`, `spectral_weights()`, `design_variance()`, `design_effect()`, `worst_case_factor()`);
- estimates variances (`random_group_variance()`, `regression_variance()`);
- evaluates samples from any design (`balance_check()`, `main_effect_share()`, `worst_case_bound()`).

The theory is in the manuscript *Exactly balanced sampling on overlapping stratifications with affine flats*
(M. Laib). The code that reproduces every result of the manuscript is in a separate repository.

Author and maintainer: Mohamed Laib, Luxembourg Institute of Science and Technology (LIST).
Contact: laib.med@gmail.com, mohamed.laib@list.lu.

```r
library(cosetsampling)
g <- grid_frame(61)
y <- as.vector(t(volcano[14:74, ]))
S <- strata_grid(g, c("rows", "columns", "antidiagonals", "diagonals"))
d <- coset_design(g, S, h = 1)
design_effect(d, y)
draw_sample(d)$units
```

Status: development version; the licence has not been chosen yet.

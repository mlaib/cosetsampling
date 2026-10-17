test_that("a flat is balanced exactly when its direction is admissible (Proposition 1)", {
  f <- coset_frame(3, 3)
  S <- strata_margins(f)
  set.seed(1)
  for (t in 1:40) {
    D <- cosetsampling:::random_subspace(3, 2, 3)
    lab <- cosetsampling:::coset_labels(f, D)
    s <- which(lab == sample(0:8, 1))
    expect_equal(balance_check(f, S, s)$exact, is_admissible(f, S, D))
  }
})

test_that("numbers of admissible directions agree with the paper (Table 2)", {
  f <- coset_frame(2, 7); S <- strata_margins(f)
  expect_equal(length(admissible_directions(f, S, 2)), 364)
  expect_equal(length(admissible_directions(f, S, 3)), 4811)
  g <- coset_frame(3, 5); T <- strata_margins(g)
  expect_equal(length(admissible_directions(g, T, 1)), 16)
  expect_equal(length(admissible_directions(g, T, 2)), 680)
})

test_that("variance, inclusion probabilities and weights agree with enumeration of the support", {
  f <- coset_frame(3, 3); S <- strata_margins(f, vars = 1:2)
  set.seed(2); y <- rnorm(f$N)
  for (r in 1:2) {
    d <- coset_design(f, S, h = 1, r = r)
    sup <- support(d)
    p <- vapply(sup, `[[`, 1, "p"); m <- vapply(sup, function(x) mean(y[x$units]), 1)
    expect_equal(sum(p), 1)
    expect_equal(sum(p * m), mean(y))
    expect_equal(sum(p * (m - mean(y))^2), design_variance(d, y))
    P <- matrix(0, f$N, f$N)
    for (x in sup) { ind <- numeric(f$N); ind[x$units] <- 1; P <- P + x$p * outer(ind, ind) }
    expect_equal(P, joint_inclusion(d, matrix = TRUE))
    w <- spectral_weights(d)
    expect_equal(sum(w), f$N / d$k - 1)
    expect_equal(sum(w[main_effect_set(f, S)]), 0)
    for (x in sup) expect_true(balance_check(f, S, x$units)$exact)
  }
})

test_that("one stratification gives the inclusion probabilities of stratified sampling (Proposition 5(b))", {
  f <- coset_frame(3, 3); S <- strata_margins(f, vars = 1)
  d <- coset_design(f, S, h = 2, r = 1)
  P <- joint_inclusion(d, matrix = TRUE)
  lab <- stratum_labels(f, S[[1]])
  nh <- d$k / 3; Nh <- f$N / 3
  same <- outer(lab, lab, "==")
  expect_equal(P[same & !diag(f$N)], rep(nh / Nh * (nh - 1) / (Nh - 1), sum(same) - f$N))
  expect_equal(P[!same], rep((nh / Nh)^2, sum(!same)))
})

test_that("equal-weight designs have design effect 1 - adjusted R^2 (Corollary 3)", {
  g <- grid_frame(7); S <- strata_grid(g, c("rows", "columns", "antidiagonals"))
  d <- coset_design(g, S, h = 1)
  expect_true(is_equal_weight(d))
  set.seed(3); y <- rnorm(g$N) + rep(1:7, each = 7)
  pw <- fourier_power(g, y); M <- main_effect_set(g, S); R <- !M; R[1] <- FALSE
  adjR2 <- 1 - (sum(pw[R]) / sum(R)) / (sum(pw[-1]) / (g$N - 1))
  expect_equal(design_effect(d, y), 1 - adjR2)
})

test_that("grid design effects agree with the paper (Table 5)", {
  y <- as.vector(t(datasets::volcano[14:74, ]))          # unit (i, j) is number 61 i + j + 1
  g <- grid_frame(61)
  fam <- c("rows", "columns", "antidiagonals", "diagonals")
  deff <- vapply(1:4, function(j) design_effect(coset_design(g, strata_grid(g, fam[1:j]), h = 1), y), 1)
  expect_equal(round(deff, 4), c(0.7171, 0.1621, 0.1197, 0.0878))
})

test_that("the random-group estimator is unbiased (Proposition 6(b))", {
  f <- coset_frame(3, 3); S <- strata_margins(f, vars = 1:2)
  d <- coset_design(f, S, h = 1, r = 3)
  set.seed(4); y <- rnorm(f$N)
  ev <- 0
  for (i in seq_along(d$directions)) {
    lab <- cosetsampling:::coset_labels(f, d$directions[[i]])
    sets <- utils::combn(d$Q, d$r, simplify = FALSE)
    for (cs in sets) {
      u <- which(lab %in% (cs - 1))
      v <- random_group_variance(list(units = u, unit_coset = lab[u]), y, d)$variance
      ev <- ev + d$probs[i] / length(sets) * v
    }
  }
  expect_equal(ev, design_variance(d, y))
})

test_that("main-effect weights from the Fourier transform agree with the stratum counts", {
  f <- coset_frame(2, 7); S <- strata_margins(f)
  set.seed(5); k <- 16
  s <- sample.int(f$N, k)
  counts <- vapply(S, function(Wp) sum(tabulate(stratum_labels(f, Wp)[s] + 1, 2)^2), 1)
  expect_equal(main_effect_weight(f, S, s, k), sum(2 * counts - k^2) / k^2)
  d <- coset_design(f, S, h = 4)
  expect_equal(main_effect_share(f, S, list(draw_sample(d)$units), k)[["tau"]], 0)
  srs <- replicate(4000, sample.int(f$N, k), simplify = FALSE)
  tau <- main_effect_share(f, S, srs, k)
  expect_lt(abs(tau[["tau"]] - 7 / 127), 4 * tau[["se"]])     # simple random sampling has tau = d_S/(N - 1)
})

test_that("the minimax mixture agrees with the paper (Table 2)", {
  skip_if_not_installed("lpSolve")
  f <- coset_frame(2, 7); S <- strata_margins(f)
  mm <- minimax_mixture(f, S, h = 2)
  expect_equal(mm$value, 5 / 18, tolerance = 1e-8)
  expect_equal(round(mm$ratio, 3), 1.075)
  expect_equal(round(worst_case_factor(mm$design), 3), 1.075)
})

test_that("samples drawn by rejection are balanced and of size k", {
  f <- coset_frame(2, 7); S <- strata_margins(f, order = 2)
  d <- coset_design(f, S, h = 4, nu = "rejection")
  set.seed(6)
  for (s in draw_sample(d, 5)) {
    expect_equal(length(s$units), 16)
    expect_true(balance_check(f, S, s$units)$exact)
  }
})

test_that("the regression estimator is unbiased for one stratification (stratified sampling)", {
  g <- grid_frame(5); S <- strata_grid(g, "rows")
  d <- coset_design(g, S, h = 1, r = 2)
  set.seed(7); y <- rnorm(g$N) + rep(1:5, each = 5)
  ev <- 0
  for (i in seq_along(d$directions)) {
    lab <- cosetsampling:::coset_labels(g, d$directions[[i]])
    sets <- utils::combn(d$Q, d$r, simplify = FALSE)
    for (cs in sets) {
      u <- which(lab %in% (cs - 1))
      ev <- ev + d$probs[i] / length(sets) * regression_variance(u, y, g, S)$variance
    }
  }
  expect_equal(ev, design_variance(d, y))
})

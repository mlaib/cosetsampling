test_that("subspace enumeration matches the Gaussian binomial", {
  for (cfg in list(c(4, 2, 3), c(5, 2, 2), c(3, 1, 5))) {
    expect_equal(length(cosetsampling:::all_subspaces(cfg[1], cfg[2], cfg[3])), gaussian_binomial(cfg[1], cfg[2], cfg[3]))
  }
})

test_that("rank and null space over GF(p)", {
  A <- matrix(c(1, 2, 3, 2, 4, 6), 2, byrow = TRUE)
  expect_equal(rank_mod(A, 7), 1)
  N <- nullspace_mod(A, 7)
  expect_equal(nrow(N), 2)
  expect_true(all((A %*% t(N)) %% 7 == 0))
})

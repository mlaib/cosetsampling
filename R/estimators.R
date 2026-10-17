#' Variance estimators for coset designs
#'
#' `random_group_variance()` is the unbiased estimator `(1 - r/Q) s_c^2 / r` of the variance of the sample mean
#' under a coset design with `r >= 2` cosets, where `s_c^2` is the sample variance of the means of the selected
#' cosets, with the Student interval on `r - 1` degrees of freedom. It is unbiased but has only `r - 1` degrees
#' of freedom.
#'
#' `regression_variance()` is the residual-based estimator `(1/k - 1/N) RSS/(k - rank)` from the least-squares
#' fit, on the sample, of the main-effects model of the stratifications, with the Student interval on
#' `k - rank` degrees of freedom. It is more stable, but not unbiased in general.
#' @param sample A sample drawn by [draw_sample()].
#' @param y The response of the whole frame, or of the selected units in the order of `sample$units`.
#' @param design The `coset_design` of the sample.
#' @param frame,strata The frame and stratifications.
#' @param level The confidence level.
#' @return A list with the estimate, its degrees of freedom and the interval for the population mean.
#' @export
random_group_variance <- function(sample, y, design, level = 0.95) {
  r <- design$r
  if (r < 2) stop("the random-group estimator needs r >= 2 cosets")
  ys <- if (length(y) == design$frame$N) y[sample$units] else y
  means <- tapply(ys, sample$unit_coset, mean)
  v <- (1 - r / design$Q) * stats::var(as.vector(means)) / r
  est <- mean(ys)
  half <- stats::qt(1 - (1 - level) / 2, r - 1) * sqrt(v)
  list(estimate = est, variance = v, df = r - 1, interval = c(est - half, est + half))
}

#' @rdname random_group_variance
#' @export
regression_variance <- function(sample, y, frame, strata, level = 0.95) {
  units <- if (is.list(sample)) sample$units else sample
  ys <- if (length(y) == frame$N) y[units] else y
  k <- length(units)
  Xs <- cbind(1, do.call(cbind, lapply(strata, function(Wp) {
    lab <- factor(stratum_labels(frame, Wp)[units], levels = 0:(frame$q^nrow(Wp) - 1))
    stats::model.matrix(~ lab)[, -1, drop = FALSE]
  })))
  qrX <- qr(Xs, tol = 1e-9)
  df <- k - qrX$rank
  if (df < 1) stop("the main-effects model has no residual degree of freedom on this sample")
  rss <- sum(qr.resid(qrX, ys)^2)
  v <- (1 / k - 1 / frame$N) * rss / df
  est <- mean(ys)
  half <- stats::qt(1 - (1 - level) / 2, df) * sqrt(v)
  list(estimate = est, variance = v, df = df, interval = c(est - half, est + half))
}

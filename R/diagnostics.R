#' Balance of a sample
#'
#' Counts the units of a sample in every stratum of every stratification and compares the counts with the
#' targets `|s| / q^n` of proportional allocation.
#' @param frame A `coset_frame`.
#' @param strata A `coset_strata` list.
#' @param units The units of the sample (frame numbers).
#' @return A list with `exact` (the sample is balanced on every stratification), `max_deviation` (the
#'   largest deviation of a count from its target, per stratification) and the `counts`.
#' @export
balance_check <- function(frame, strata, units) {
  k <- length(units)
  res <- lapply(strata, function(Wp) {
    nlev <- frame$q^nrow(Wp)
    cnt <- tabulate(stratum_labels(frame, Wp)[units] + 1, nlev)
    list(counts = cnt, dev = max(abs(cnt - k / nlev)))
  })
  dev <- vapply(res, `[[`, 1, "dev")
  list(exact = all(dev < 1e-9), max_deviation = dev, counts = lapply(res, `[[`, "counts"))
}

#' Main-effect weight of samples
#'
#' The main-effect weight of a sample is
#' \eqn{t_M(s) = \sum_{u \in M_S} |\eta_s(u)|^2}{t_M(s) = sum_{u in M_S} |eta_s(u)|^2}, with
#' \eqn{\eta_s(u) = k^{-1} \sum_{x \in s} \exp(2 \pi i \, u \cdot x / q)}{eta_s(u) = k^-1 sum_{x in s} exp(2 pi i u.x/q)}
#' and `k` the target sample size; it is zero exactly when the
#' sample is balanced on every stratification. Its mean over the samples of a design, divided by `N/k - 1`, is
#' the main-effect weight share `tau` of the design. Samples may come from any design, for instance the cube
#' method.
#'
#' `worst_case_bound()` estimates the lower bound on the worst-case variance over the responses whose main
#' effects carry at most `C` times the residual variance, relative to an equal-weight design:
#' `1 + tau (C |R_S|/d_S - 1) - var(|s|)/{k^2 (N/k - 1)}`, where the last term, for designs whose size varies,
#' restricts the bound to responses with mean zero. Both are estimated by means over the samples, with Monte
#' Carlo standard errors.
#' @param frame A `coset_frame`.
#' @param strata A `coset_strata` list.
#' @param samples A list of samples (vectors of frame numbers).
#' @param k The target sample size, for which the inclusion probabilities are `k/N`.
#' @param C The bound on the main effects relative to the residual variance.
#' @return `main_effect_weight()`: the vector of `t_M(s)`. `main_effect_share()`: `tau` and its standard
#'   error. `worst_case_bound()`: the bound and its standard error.
#' @export
main_effect_weight <- function(frame, strata, samples, k) {
  if (!is.list(samples)) samples <- list(samples)
  M <- main_effect_set(frame, strata)
  dims <- if (frame$m == 1) NULL else rep(frame$q, frame$m)
  vapply(samples, function(s) {
    ind <- numeric(frame$N); ind[s] <- 1
    F <- stats::fft(if (is.null(dims)) ind else array(ind, dims))
    sum(Mod(as.vector(F)[M])^2) / k^2
  }, 1)
}

#' @rdname main_effect_weight
#' @export
main_effect_share <- function(frame, strata, samples, k) {
  T0 <- frame$N / k - 1
  tm <- main_effect_weight(frame, strata, samples, k)
  c(tau = mean(tm) / T0, se = stats::sd(tm) / sqrt(length(tm)) / T0)
}

#' @rdname main_effect_weight
#' @export
worst_case_bound <- function(frame, strata, samples, k, C = 1) {
  if (!is.list(samples)) samples <- list(samples)
  T0 <- frame$N / k - 1
  dr <- main_effect_count(frame, strata)
  tm <- main_effect_weight(frame, strata, samples, k)
  size <- lengths(samples)
  z <- 1 + (C * dr[["residual"]] / dr[["d"]] - 1) * tm / T0 - (size - k)^2 / (k^2 * T0)
  c(bound = mean(z), se = stats::sd(z) / sqrt(length(z)))
}

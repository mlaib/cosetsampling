#' Coset designs
#'
#' The coset design `C(h, nu, r)` draws a direction `V` of dimension `h` from the mixing distribution `nu` over
#' admissible directions, selects `r` of the `Q = q^(m - h)` cosets of `V` by simple random sampling without
#' replacement, and returns their union. Every sample has `k = r q^h` units and is exactly balanced on every
#' stratification for which the directions are admissible.
#'
#' @param frame A `coset_frame`.
#' @param strata A `coset_strata` list.
#' @param h The dimension of the flats.
#' @param r The number of cosets, `1 <= r < q^(m - h)`.
#' @param nu The mixing distribution: `"uniform"` over all admissible directions (listed by
#'   [admissible_directions()]); `"rejection"`, uniform over admissible directions drawn by rejection without
#'   listing them (only sampling is then available); or a list of dual directions, with `probs`.
#' @param probs Probabilities of the directions in `nu`, by default equal.
#' @return An object of class `coset_design`.
#' @references Laib, M. Exactly balanced sampling on overlapping stratifications with affine flats.
#'   Manuscript.
#' @examples
#' g <- grid_frame(7)
#' S <- strata_grid(g, c("rows", "columns", "antidiagonals"))
#' d <- coset_design(g, S, h = 1)      # equal-weight design: uniform over the 5 other line directions
#' d
#' s <- draw_sample(d)
#' balance_check(g, S, s$units)$exact
#' @export
coset_design <- function(frame, strata, h, r = 1, nu = "uniform", probs = NULL) {
  Q <- frame$q^(frame$m - h)
  if (r < 1 || r >= Q) stop("r must satisfy 1 <= r < q^(m - h)")
  dirs <- NULL
  if (identical(nu, "uniform")) {
    dirs <- admissible_directions(frame, strata, h)
    if (!length(dirs)) stop("no admissible direction of this dimension")
    probs <- rep(1 / length(dirs), length(dirs))
  } else if (identical(nu, "rejection")) {
    if (!isTRUE(admissible_exists(frame, strata, h))) message("admissible directions may not exist")
  } else {
    dirs <- lapply(nu, function(D) rref_mod(matrix(D, ncol = frame$m), frame$q)$R)
    if (any(vapply(dirs, nrow, 1) != frame$m - h)) stop("each direction needs m - h independent rows")
    if (!all(vapply(dirs, function(D) is_admissible(frame, strata, D), TRUE))) stop("a direction is not admissible")
    probs <- if (is.null(probs)) rep(1 / length(dirs), length(dirs)) else probs / sum(probs)
  }
  structure(list(frame = frame, strata = strata, h = h, r = r, Q = Q, k = r * frame$q^h,
                 directions = dirs, probs = probs, implicit = is.null(dirs)), class = "coset_design")
}

#' @export
print.coset_design <- function(x, ...) {
  f <- x$frame
  cat(sprintf("Coset design C(h = %d, nu, r = %d) on GF(%d)^%d: N = %d, sample size k = %d, %d stratifications\n",
              x$h, x$r, f$q, f$m, f$N, x$k, length(x$strata)))
  if (x$implicit) cat("nu: uniform over admissible directions, drawn by rejection\n")
  else cat(sprintf("nu: %d direction(s)\n", length(x$directions)))
  invisible(x)
}

need_explicit <- function(design) if (design$implicit) stop("this needs the directions of nu listed (not nu = \"rejection\")")

# label of the coset of V = D^perp containing each unit (0, ..., q^(m-h) - 1)
coset_labels <- function(frame, D) drop(((frame$X %*% t(D)) %% frame$q) %*% frame$q^((nrow(D) - 1):0))

#' Draw samples from a coset design
#'
#' @param design A `coset_design`.
#' @param n The number of samples.
#' @return For `n = 1`, a list with the selected `units` (frame numbers), the dual direction `D`, the selected
#'   `cosets` and the coset of each selected unit (`unit_coset`); for `n > 1`, a list of such lists.
#' @export
draw_sample <- function(design, n = 1) {
  one <- function() {
    f <- design$frame
    D <- if (design$implicit) sample_direction(f, design$strata, design$h)
         else design$directions[[sample.int(length(design$directions), 1, prob = design$probs)]]
    lab <- coset_labels(f, D)
    cos <- sample.int(design$Q, design$r) - 1
    units <- which(lab %in% cos)
    list(units = units, D = D, cosets = cos, unit_coset = lab[units])
  }
  if (n == 1) one() else replicate(n, one(), simplify = FALSE)
}

#' Inclusion probabilities
#'
#' Every coset design has first-order inclusion probabilities `k/N`. The second-order inclusion probability
#' of two units depends only on their difference `z`: `pi(z) = (k/N) {rho(z) + (1 - rho(z)) (r - 1)/(Q - 1)}`,
#' where `rho(z)` is the probability that `z` lies in the random direction.
#' @param design A `coset_design` with listed directions.
#' @param matrix If `TRUE`, the `N x N` matrix of joint inclusion probabilities; otherwise `pi(z)` for every
#'   difference `z`, in unit order (the first entry, `z = 0`, is `k/N`).
#' @export
joint_inclusion <- function(design, matrix = FALSE) {
  need_explicit(design)
  f <- design$frame
  rho <- numeric(f$N)
  for (i in seq_along(design$directions)) {
    inV <- rowSums((f$X %*% t(design$directions[[i]])) %% f$q) == 0
    rho[inV] <- rho[inV] + design$probs[i]
  }
  pz <- design$k / f$N * (rho + (1 - rho) * (design$r - 1) / (design$Q - 1))
  pz[1] <- design$k / f$N
  if (!matrix) return(pz)
  if (f$N > 6000) stop("N is too large for the full matrix")
  code <- matrix(0, f$N, f$N)
  for (j in seq_len(f$m)) code <- code + (outer(f$X[, j], f$X[, j], "-") %% f$q) * f$q^(f$m - j)
  matrix(pz[code + 1], f$N, f$N)
}

#' @rdname joint_inclusion
#' @export
inclusion_probabilities <- function(design) rep(design$k / design$frame$N, design$frame$N)

#' Spectral weights
#'
#' The variance of the Horvitz-Thompson estimator of the mean under a translation-invariant design is
#' `sum_u w(u) |yhat(u)|^2` over the nonzero frequencies `u`. For a coset design
#' `w(u) = c(u) (Q - r)/(r (Q - 1))`, with `c(u)` the probability that `u` lies in the random dual direction.
#' The weights sum to `N/k - 1`; an exactly balanced design puts none on the main-effect frequencies.
#' @param design A `coset_design` with listed directions.
#' @return The weights `w(u)` in unit order, with `w(0) = 0`.
#' @export
spectral_weights <- function(design) {
  need_explicit(design)
  f <- design$frame
  cov <- numeric(f$N)
  for (i in seq_along(design$directions)) {
    sp <- span_codes(design$directions[[i]], f$q, f$m) + 1
    cov[sp] <- cov[sp] + design$probs[i]
  }
  w <- cov * (design$Q - design$r) / (design$r * (design$Q - 1))
  w[1] <- 0
  w
}

#' Fourier power of a response
#'
#' `|yhat(u)|^2` for every frequency, where `yhat(u) = N^-1 sum_x y(x) exp(-2 pi i u.x / q)`.
#' @param frame A `coset_frame`.
#' @param y A response, one value per unit in frame order.
#' @export
fourier_power <- function(frame, y) {
  if (length(y) != frame$N) stop("y needs one value per unit")
  Y <- if (frame$m == 1) y else array(y, rep(frame$q, frame$m))
  as.vector(Mod(stats::fft(Y) / frame$N)^2)
}

#' Exact variance and design effect
#'
#' `design_variance()` is the exact variance of the Horvitz-Thompson estimator of the mean (the sample mean)
#' under a coset design; `design_effect()` divides it by the variance under simple random sampling of the same
#' size; `srs_variance()` is that variance.
#' @param design A `coset_design` with listed directions.
#' @param y A response, one value per unit in frame order.
#' @export
design_variance <- function(design, y) sum(spectral_weights(design) * fourier_power(design$frame, y))

#' @rdname design_variance
#' @param N,k Population and sample sizes.
#' @export
srs_variance <- function(y, N = length(y), k) (1 / k - 1 / N) * stats::var(y)

#' @rdname design_variance
#' @export
design_effect <- function(design, y) design_variance(design, y) / srs_variance(y, design$frame$N, design$k)

#' Worst-case factor
#'
#' The worst-case variance of an exactly balanced translation-invariant design over the responses whose
#' residual variance from the main-effects model is at most `sigma^2` is `sigma^2 max_{R_S} w(u)`. The
#' worst-case factor divides it by the lower bound `(N/k - 1)/|R_S|`, which equal-weight designs attain; it is
#' 1 for them.
#' @param design A `coset_design` with listed directions.
#' @param tol Tolerance for `is_equal_weight()`.
#' @export
worst_case_factor <- function(design) {
  f <- design$frame
  R <- !main_effect_set(f, design$strata); R[1] <- FALSE
  max(spectral_weights(design)[R]) / ((f$N / design$k - 1) / sum(R))
}

#' @rdname worst_case_factor
#' @export
is_equal_weight <- function(design, tol = 1e-10) abs(worst_case_factor(design) - 1) < tol

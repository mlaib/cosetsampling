#' Admissible directions
#'
#' A flat `a + V` is balanced on the stratification `W` exactly when `V + W` is the whole space, that is when
#' `V^perp` and `W^perp` meet only in 0. A direction is admissible when it is transversal to every
#' stratification. Directions are stored as a basis of the dual direction `D = V^perp`, a subspace of
#' dimension `m - h` for flats of dimension `h` (size `q^h`).
#'
#' `admissible_directions()` lists every admissible direction by enumerating the subspaces of dimension
#' `m - h`; their number grows like \eqn{q^{h(m - h)}}{q^(h (m - h))}, so this is for moderate frames. `sample_direction()`
#' draws a uniform admissible direction by rejection, without listing them. `is_admissible()` checks one
#' direction.
#'
#' @param frame A `coset_frame`.
#' @param strata A `coset_strata` list.
#' @param h The dimension of the flats.
#' @param D A dual direction: a matrix with `m - h` rows whose row space is `V^perp`.
#' @param max_subspaces Largest number of subspaces to enumerate.
#' @param max_tries Largest number of rejection steps.
#' @return `admissible_directions()` a list of dual directions; `sample_direction()` one dual direction;
#'   `is_admissible()` a logical.
#' @examples
#' f <- coset_frame(2, 7)
#' S <- strata_margins(f)
#' length(admissible_directions(f, S, h = 2))   # 364, flats of 4 units
#' @export
admissible_directions <- function(frame, strata, h, max_subspaces = 2e5) {
  k <- frame$m - h
  if (k < 1 || h < 1) stop("h must be between 1 and m - 1")
  n <- gaussian_binomial(frame$m, k, frame$q)
  if (n > max_subspaces) stop(sprintf("%.0f subspaces to enumerate; use sample_direction() or raise max_subspaces", n))
  M <- main_effect_set(frame, strata)
  Filter(function(D) !any(M[span_codes(D, frame$q, frame$m) + 1]), all_subspaces(frame$m, k, frame$q))
}

#' @rdname admissible_directions
#' @export
is_admissible <- function(frame, strata, D) {
  M <- main_effect_set(frame, strata)
  !any(M[span_codes(D, frame$q, frame$m) + 1])
}

#' @rdname admissible_directions
#' @export
sample_direction <- function(frame, strata, h, max_tries = 1e5) {
  M <- main_effect_set(frame, strata)
  for (t in seq_len(max_tries)) {
    D <- random_subspace(frame$m, frame$m - h, frame$q)
    if (!any(M[span_codes(D, frame$q, frame$m) + 1])) return(D)
  }
  stop("no admissible direction found; there may be none (see admissible_exists())")
}

#' Existence of admissible directions
#'
#' Checks the sufficient condition of Theorem 1(b) of the accompanying paper: when the number of projective
#' points of `M_S` is below the number of points of a projective space of dimension `h`, an admissible
#' direction exists. Otherwise the answer is found by enumeration when it is feasible.
#' @inheritParams admissible_directions
#' @return `TRUE` or `FALSE`, with an attribute `reason`.
#' @export
admissible_exists <- function(frame, strata, h, max_subspaces = 2e5) {
  q <- frame$q
  npts <- sum(main_effect_set(frame, strata)) / (q - 1)
  theta <- (q^(h + 1) - 1) / (q - 1)
  if (npts < theta) return(structure(TRUE, reason = "sufficient condition of Theorem 1(b)"))
  if (gaussian_binomial(frame$m, frame$m - h, q) <= max_subspaces)
    return(structure(length(admissible_directions(frame, strata, h, max_subspaces)) > 0, reason = "enumeration"))
  structure(NA, reason = "undecided: too many subspaces to enumerate")
}

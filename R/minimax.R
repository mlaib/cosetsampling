#' Minimax mixture of directions
#'
#' Solves the linear program that minimises `max_{u in R_S} c(u)` over probability distributions `nu` on the
#' given admissible directions, where `c(u)` is the probability that the residual frequency `u` lies in the
#' random dual direction. The result is the minimax design within the class of coset designs `C(h, nu, r)` for
#' fixed `h` and `r` (the program does not depend on `r`). Its value is at least `(Q - 1)/|R_S|`, with
#' equality exactly when an equal-weight coset design with flats of dimension `h` exists. Needs the package
#' lpSolve.
#' @param frame A `coset_frame`.
#' @param strata A `coset_strata` list.
#' @param h The dimension of the flats.
#' @param directions The candidate directions, by default all admissible ones.
#' @param r The number of cosets of the returned design.
#' @return A list with the `design` (a `coset_design` on the directions with positive probability), the
#'   program's `value`, the lower `bound` `(Q - 1)/|R_S|` and their `ratio`.
#' @examples
#' if (requireNamespace("lpSolve", quietly = TRUE)) {
#'   f <- coset_frame(2, 5)
#'   mm <- minimax_mixture(f, strata_margins(f), h = 2)
#'   mm$ratio
#' }
#' @export
minimax_mixture <- function(frame, strata, h, directions = admissible_directions(frame, strata, h), r = 1) {
  if (!requireNamespace("lpSolve", quietly = TRUE)) stop("minimax_mixture() needs the package lpSolve")
  if (!length(directions)) stop("no admissible direction")
  R <- !main_effect_set(frame, strata); R[1] <- FALSE
  rows <- which(R)
  A <- matrix(0, length(rows), length(directions))
  for (j in seq_along(directions)) {
    sp <- span_codes(directions[[j]], frame$q, frame$m) + 1
    A[, j] <- rows %in% sp
  }
  nd <- length(directions)
  sol <- lpSolve::lp("min", c(rep(0, nd), 1),
                     rbind(cbind(A, -1), c(rep(1, nd), 0)),
                     c(rep("<=", nrow(A)), "="), c(rep(0, nrow(A)), 1))
  if (sol$status != 0) stop("the linear program was not solved")
  nu <- sol$solution[seq_len(nd)]
  keep <- nu > 1e-12
  Q <- frame$q^(frame$m - h)
  bound <- (Q - 1) / length(rows)
  list(design = coset_design(frame, strata, h, r, nu = directions[keep], probs = nu[keep]),
       value = sol$objval, bound = bound, ratio = sol$objval / bound)
}

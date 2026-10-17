#' Stratifications by cosets of subspaces
#'
#' A stratification is given by a subspace `W` of `GF(q)^m`; its strata are the cosets `a + W`. It is stored
#' as a basis of the annihilator `W^perp`, the linear forms that are constant on each stratum: the stratum of
#' unit `x` is the value of `W^perp x`.
#'
#' `strata_linear()` takes the linear forms directly, one stratification per matrix. `strata_margins()`
#' gives the one-variable stratifications (`order = 1`) or, in addition, all two-variable cross-classifications
#' (`order = 2`) of a factorial frame. `strata_grid()` gives line stratifications of a square grid: rows
#' (`i` constant), columns (`j` constant), anti-diagonals (`i + j` constant) and diagonals (`i - j` constant),
#' all wrapped around the grid, or the lines `j - c i` constant for given slopes `c`.
#'
#' @param frame A `coset_frame`.
#' @param forms A list of matrices with `m` columns, or one matrix per stratification.
#' @param order 1 for the one-variable margins, 2 to add the two-variable cross-classifications.
#' @param vars The variables to use.
#' @param which Line families of the grid.
#' @param slopes Further slopes `c` of lines `j - c i = b`.
#' @return A list of annihilator bases, of class `coset_strata`.
#' @examples
#' g <- grid_frame(5)
#' S <- strata_grid(g, c("rows", "columns"))
#' main_effect_count(g, S)
#' @export
strata_linear <- function(frame, forms) {
  if (is.matrix(forms)) forms <- list(forms)
  out <- lapply(forms, function(L) {
    L <- matrix(L, ncol = frame$m) %% frame$q
    e <- rref_mod(L, frame$q)
    if (!length(e$pivots)) stop("a stratification needs at least one nonzero linear form")
    e$R
  })
  structure(out, class = "coset_strata")
}

#' @rdname strata_linear
#' @export
strata_margins <- function(frame, order = 1, vars = seq_len(frame$m)) {
  E <- diag(frame$m)
  forms <- lapply(vars, function(i) E[i, , drop = FALSE])
  if (order >= 2) forms <- c(forms, lapply(utils::combn(vars, 2, simplify = FALSE), function(ij) E[ij, , drop = FALSE]))
  strata_linear(frame, forms)
}

#' @rdname strata_linear
#' @export
strata_grid <- function(frame, which = c("rows", "columns", "antidiagonals", "diagonals"), slopes = NULL) {
  if (frame$m != 2) stop("strata_grid() needs a grid frame (m = 2)")
  q <- frame$q
  dict <- list(rows = c(1, 0), columns = c(0, 1), antidiagonals = c(1, 1), diagonals = c(1, q - 1))
  which <- match.arg(which, names(dict), several.ok = TRUE)
  forms <- lapply(dict[which], function(v) matrix(v, 1))
  if (length(slopes)) forms <- c(forms, lapply(slopes, function(cc) matrix(c(-cc, 1) %% q, 1)))
  strata_linear(frame, forms)
}

#' Stratum labels
#'
#' @param frame A `coset_frame`.
#' @param Wp The annihilator basis of one stratification.
#' @return For each unit, the label (`0, ..., q^n - 1`) of its stratum.
#' @export
stratum_labels <- function(frame, Wp) {
  Wp <- matrix(Wp, ncol = frame$m)
  drop(((frame$X %*% t(Wp)) %% frame$q) %*% frame$q^((nrow(Wp) - 1):0))
}

#' Main-effect and residual frequencies
#'
#' The main-effect frequencies `M_S` are the nonzero vectors of the union of the annihilators `W^perp`; they
#' carry the variance explained by the main-effects model of the stratifications. The residual frequencies
#' `R_S` are the other nonzero vectors.
#' @param frame A `coset_frame`.
#' @param strata A `coset_strata` list.
#' @return `main_effect_set()` a logical vector over the frequencies (in unit order); `main_effect_count()` the
#'   numbers `d_S = |M_S|` and `|R_S|`.
#' @export
main_effect_set <- function(frame, strata) {
  M <- logical(frame$N)
  for (Wp in strata) M[span_codes(Wp, frame$q, frame$m) + 1] <- TRUE
  M[1] <- FALSE
  M
}

#' @rdname main_effect_set
#' @export
main_effect_count <- function(frame, strata) {
  d <- sum(main_effect_set(frame, strata))
  c(d = d, residual = frame$N - 1 - d)
}

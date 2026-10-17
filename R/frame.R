#' Frames indexed by a vector space over a prime field
#'
#' A coset frame has `N = q^m` units, one for each vector `x` of `GF(q)^m`, with `q` prime. Units are
#' numbered `1, ..., N` in the order of their code `sum_j x_j q^(m - j)`, so that the last coordinate varies
#' fastest.
#'
#' @param q A prime: the number of categories of each variable, or the side of a square grid.
#' @param m The dimension (number of variables).
#' @return An object of class `coset_frame` with elements `q`, `m`, `N` and `X`, the `N x m` matrix of
#'   coordinates in `0, ..., q - 1`.
#' @examples
#' f <- coset_frame(3, 2)
#' f$X
#' @export
coset_frame <- function(q, m) {
  q <- as.integer(q); m <- as.integer(m)
  if (q < 2 || any(q %% seq_len(max(1, floor(sqrt(q))))[-1] == 0)) stop("q must be a prime")
  if (m < 1) stop("m must be at least 1")
  X <- digits_matrix(m, q)
  structure(list(q = q, m = m, N = q^m, X = X), class = "coset_frame")
}

#' @export
print.coset_frame <- function(x, ...) {
  cat(sprintf("Coset frame: GF(%d)^%d, N = %d units\n", x$q, x$m, x$N))
  invisible(x)
}

#' Square grid of prime side
#'
#' The `q x q` grid as the frame `GF(q)^2`: unit `(i, j)`, row `i` and column `j` in `0, ..., q - 1`, is unit
#' number `q * i + j + 1`.
#' @param q A prime.
#' @return A `coset_frame` with `m = 2`.
#' @export
grid_frame <- function(q) coset_frame(q, 2)

#' Complete factorial frame from data
#'
#' Maps a complete factorial population, one unit for every combination of the categories of `m`
#' variables with `q` categories each (`q` prime), to the frame `GF(q)^m`.
#' @param data A data frame whose columns are the `m` variables (factors or vectors), with `q^m` rows.
#' @return A list with the `frame` and `unit`, the frame unit of each row of `data`.
#' @export
factorial_frame <- function(data) {
  data <- as.data.frame(data)
  lev <- lapply(data, function(v) sort(unique(as.character(v))))
  q <- unique(lengths(lev))
  if (length(q) != 1) stop("all variables must have the same number of categories")
  m <- ncol(data)
  fr <- coset_frame(q, m)
  if (nrow(data) != fr$N) stop("a complete factorial frame has q^m rows")
  codes <- mapply(function(v, l) match(as.character(v), l) - 1, data, lev)
  codes <- matrix(codes, ncol = m)
  unit <- drop(codes %*% q^((m - 1):0)) + 1
  if (anyDuplicated(unit)) stop("every combination of categories must occur exactly once")
  list(frame = fr, unit = unit, levels = lev)
}

unit_codes <- function(frame) seq_len(frame$N) - 1

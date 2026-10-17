# Linear algebra over the prime field GF(p). Matrices hold integers in 0..p-1 stored as doubles, which is exact
# for the moduli used here (p^2 well below 2^53).

inv_mod <- function(a, p) {
  # inverse of a (nonzero mod p) by the extended Euclidean algorithm
  a <- a %% p
  t0 <- 0; t1 <- 1; r0 <- p; r1 <- a
  while (r1 != 0) {
    qt <- r0 %/% r1
    tmp <- r0 - qt * r1; r0 <- r1; r1 <- tmp
    tmp <- t0 - qt * t1; t0 <- t1; t1 <- tmp
  }
  if (r0 != 1) stop("not invertible modulo p")
  t0 %% p
}

#' Reduced row echelon form over a prime field
#'
#' @param A Integer matrix.
#' @param p A prime.
#' @return A list with the nonzero rows `R` of the reduced row echelon form of `A` modulo `p` and the
#'   pivot columns `pivots`.
#' @examples
#' rref_mod(matrix(c(1, 2, 2, 4), 2), 5)
#' @export
rref_mod <- function(A, p) {
  R <- as.matrix(A) %% p
  storage.mode(R) <- "double"
  nr <- nrow(R); nc <- ncol(R)
  pivots <- integer(0); r <- 1L
  for (cc in seq_len(nc)) {
    if (r > nr) break
    nz <- which(R[r:nr, cc] != 0)
    if (!length(nz)) next
    i <- r + nz[1] - 1L
    if (i != r) { tmp <- R[r, ]; R[r, ] <- R[i, ]; R[i, ] <- tmp }
    R[r, ] <- (R[r, ] * inv_mod(R[r, cc], p)) %% p
    for (i2 in seq_len(nr)) if (i2 != r && R[i2, cc] != 0) R[i2, ] <- (R[i2, ] - R[i2, cc] * R[r, ]) %% p
    pivots <- c(pivots, cc); r <- r + 1L
  }
  list(R = R[seq_len(r - 1L), , drop = FALSE], pivots = pivots)
}

#' Rank over a prime field
#' @inheritParams rref_mod
#' @return The rank of `A` modulo `p`.
#' @export
rank_mod <- function(A, p) if (length(A)) length(rref_mod(A, p)$pivots) else 0L

#' Null space over a prime field
#' @inheritParams rref_mod
#' @return A matrix whose rows form a basis of `{x : A x = 0 (mod p)}`.
#' @export
nullspace_mod <- function(A, p) {
  A <- as.matrix(A)
  m <- ncol(A)
  e <- rref_mod(A, p)
  free <- setdiff(seq_len(m), e$pivots)
  out <- matrix(0, length(free), m)
  for (f in seq_along(free)) {
    out[f, free[f]] <- 1
    for (i in seq_along(e$pivots)) out[f, e$pivots[i]] <- (-e$R[i, free[f]]) %% p
  }
  out
}

# all coefficient vectors of length d over GF(p), as rows (p^d x d), first coordinate most significant
digits_matrix <- function(d, p) {
  n <- p^d
  if (d == 0) return(matrix(0, 1, 0))
  idx <- 0:(n - 1)
  out <- matrix(0, n, d)
  for (j in d:1) { out[, j] <- idx %% p; idx <- idx %/% p }
  out
}

# codes (0-based) of the subspace spanned by the rows of B in GF(p)^m
span_codes <- function(B, p, m) {
  B <- matrix(B, ncol = m)
  if (nrow(B) == 0) return(0)
  V <- (digits_matrix(nrow(B), p) %*% B) %% p
  sort(unique(drop(V %*% p^((m - 1):0))))
}

#' Number of subspaces of a given dimension
#'
#' The Gaussian binomial coefficient: the number of `k`-dimensional subspaces of `GF(q)^m`.
#' @param m,k Dimensions.
#' @param q Field size.
#' @export
gaussian_binomial <- function(m, k, q) {
  if (k < 0 || k > m) return(0)
  num <- 1; den <- 1
  for (i in seq_len(k) - 1) { num <- num * (q^(m - i) - 1); den <- den * (q^(i + 1) - 1) }
  num / den
}

# reduced row echelon bases (k x m) of all k-dimensional subspaces of GF(p)^m
all_subspaces <- function(m, k, p) {
  if (k == 0) return(list(matrix(0, 0, m)))
  out <- list()
  for (piv in utils::combn(m, k, simplify = FALSE)) {
    free <- do.call(rbind, lapply(seq_len(k), function(i) {
      js <- setdiff(seq.int(piv[i], m), piv)
      js <- js[js > piv[i]]
      if (length(js)) cbind(i, js) else NULL
    }))
    nfree <- if (is.null(free)) 0 else nrow(free)
    vals <- digits_matrix(nfree, p)
    for (v in seq_len(nrow(vals))) {
      B <- matrix(0, k, m)
      B[cbind(seq_len(k), piv)] <- 1
      if (nfree) B[free] <- vals[v, ]
      out[[length(out) + 1L]] <- B
    }
  }
  out
}

# a uniformly random k-dimensional subspace of GF(p)^m, as a reduced row echelon basis
random_subspace <- function(m, k, p) {
  repeat {
    A <- matrix(sample.int(p, k * m, replace = TRUE) - 1, k, m)
    e <- rref_mod(A, p)
    if (length(e$pivots) == k) return(e$R)
  }
}

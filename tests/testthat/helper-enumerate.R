# Enumerate the support of a coset design with listed directions: list of (probability, units)
support <- function(design) {
  f <- design$frame
  out <- list()
  for (i in seq_along(design$directions)) {
    lab <- cosetsampling:::coset_labels(f, design$directions[[i]])
    sets <- utils::combn(design$Q, design$r, simplify = FALSE)
    for (cs in sets) out[[length(out) + 1]] <- list(p = design$probs[i] / length(sets), units = which(lab %in% (cs - 1)))
  }
  out
}

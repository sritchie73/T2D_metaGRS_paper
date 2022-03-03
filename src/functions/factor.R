# base::factor, but allows you to explicitly set the reference group (first factor level)
# without having to list out all the levels
factor <- function(..., reference) {
  # Make sure ... arguments are named as in base::factor so we can extract/check relevant
  # args
  posargnames <- c("x", "levels", "labels", "exclude", "ordered", "nmax")
  arglist <- list(...)
  if (is.null(names(arglist))) {
    names(arglist) <- posargnames[seq_along(arglist)]
  } else {
    unnamed <- which(names(arglist) == "")
    arglist[unnamed] <- posargnames[unnamed]
  }

  if (missing(reference)) {
    base::factor(...)
  } else {
    stopifnot(length(reference) == 1)
    stopifnot(!("levels" %in% names(arglist)))
    base::factor(..., levels=c(reference, setdiff(unique(arglist$x), reference)))
  }
}

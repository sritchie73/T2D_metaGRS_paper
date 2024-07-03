library(Rmpfr)

# Two-sided
pval_from_zscore <- function(Z, precBits=53, format=TRUE) { # precBits=53 equivalent to precision of a double
  extreme_p <- 2*Rmpfr::pnorm(mpfr(Z, precBits=precBits), lower.tail=FALSE, log.p=FALSE)
  if (format) {
    extreme_p <- sapply(extreme_p, function(this_p) {
      if (round(this_p, digits=3) > 0.009) {
        suppressWarnings(format(this_p, digits=2, decimal.mark=""))
      } else {
        suppressWarnings(format(this_p, digits=1, decimal.mark=""))
      }
    })
  }
  return(extreme_p)
}

# Based on log HR and SE
pval_from_logHR_and_SE <- function(logHR, SE, ...) {
  pval_from_zscore(logHR/SE, ...)
}



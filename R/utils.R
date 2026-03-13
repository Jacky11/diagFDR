#' @importFrom utils head tail
NULL

#' @importFrom rlang .data
NULL

#' @keywords internal
`%||%` <- function(a, b) {
  if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
}

#' @keywords internal
safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

#' @keywords internal
jaccard_vec <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  if (length(a) == 0 || length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

#' @keywords internal
wilson_ci <- function(x, n, conf = 0.95) {
  if (n == 0) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - (1 - conf) / 2)
  phat <- x / n
  denom <- 1 + z^2 / n
  center <- (phat + z^2 / (2*n)) / denom
  half <- (z * sqrt(phat*(1 - phat)/n + z^2/(4*n^2))) / denom
  c(center - half, center + half)
}

#' @keywords internal
choose_low_conf <- function(qmax,
                            low_default = c(0.2, 0.5),
                            low_fallback = c(0.1, 0.2)) {
  if (is.finite(qmax) && qmax >= low_default[2] - 1e-12) return(low_default)
  c(low_fallback[1], min(low_fallback[2], qmax))
}

#' @keywords internal
dfdr_safe_filename <- function(x) {
  # replace anything not safe for filenames with "_"
  gsub("[^A-Za-z0-9._-]+", "_", x)
}


#' Safe minimum that handles NA values
#'
#' Returns the minimum of x, with NA values removed.
#' If all values are NA, returns NA_real_.
#'
#' @param x Numeric vector
#' @return Minimum value or NA_real_
#' @keywords internal
safe_min <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

#' Safe maximum that handles NA values
#'
#' Returns the maximum of x, with NA values removed.
#' If all values are NA, returns NA_real_.
#'
#' @param x Numeric vector
#' @return Maximum value or NA_real_
#' @keywords internal
safe_max <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

#' Safe mean that handles NA values
#'
#' Returns the mean of x, with NA values removed.
#' If all values are NA, returns NA_real_.
#'
#' @param x Numeric vector
#' @return Mean value or NA_real_
#' @keywords internal
safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}


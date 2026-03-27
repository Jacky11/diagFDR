#' @importFrom utils head tail
NULL

#' @importFrom rlang .data
NULL

#' Internal: "or else" helper for defaults
#'
#' Returns \code{a} if it is non-\code{NULL} and contains at least one non-\code{NA}
#' value; otherwise returns \code{b}. This is used throughout the package to
#' provide fallbacks for optional metadata fields.
#'
#' @param a,b Objects to choose from.
#'
#' @return \code{a} if available, otherwise \code{b}.
#'
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) {
  if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
}

#' Internal: safe finite minimum
#'
#' Returns the minimum of finite values in \code{x}. Non-finite values
#' (\code{NA}, \code{Inf}, \code{-Inf}) are removed. If no finite values remain,
#' returns \code{NA_real_}.
#'
#' @param x Numeric vector.
#'
#' @return Numeric scalar minimum of finite values, or \code{NA_real_}.
#'
#' @keywords internal
safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

#' Internal: Jaccard similarity for two ID vectors
#'
#' Computes the Jaccard index \eqn{|A \cap B|/|A \cup B|} on unique values of \code{a}
#' and \code{b}. Returns \code{NA_real_} if either set is empty (or both are empty),
#' which is convenient for diagnostics where the overlap is undefined without
#' discoveries in both sets.
#'
#' @param a,b Vectors representing sets.
#'
#' @return Numeric scalar in \eqn{[0,1]} or \code{NA_real_} if undefined.
#'
#' @keywords internal
jaccard_vec <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  if (length(a) == 0 || length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

#' Internal: Wilson score confidence interval for a binomial proportion
#'
#' Computes a Wilson score interval for \eqn{x/n} at confidence level \code{conf}.
#' Returns \code{c(NA, NA)} when \code{n == 0}.
#'
#' @param x Integer number of successes.
#' @param n Integer number of trials.
#' @param conf Confidence level in \eqn{(0,1)} (default 0.95).
#'
#' @return Numeric vector of length 2: lower and upper confidence bounds.
#'
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

#' Internal: choose a low-confidence q-interval based on export range
#'
#' Chooses an interval for pooled "low-confidence" diagnostics depending on the
#' available maximum q-value \code{qmax}. If \code{qmax} covers the default interval,
#' returns \code{low_default}; otherwise falls back to an interval capped by \code{qmax}.
#'
#' @param qmax Numeric scalar. Maximum available q-value.
#' @param low_default Length-2 numeric vector. Default low-confidence interval.
#' @param low_fallback Length-2 numeric vector. Fallback interval when \code{qmax} is smaller.
#'
#' @return Length-2 numeric vector \code{c(lo, hi)}.
#'
#' @keywords internal
choose_low_conf <- function(qmax,
                            low_default = c(0.2, 0.5),
                            low_fallback = c(0.1, 0.2)) {
  if (is.finite(qmax) && qmax >= low_default[2] - 1e-12) return(low_default)
  c(low_fallback[1], min(low_fallback[2], qmax))
}

#' Internal: make a string safe for filenames
#'
#' Replaces characters that are typically unsafe in filenames with underscores.
#' Keeps letters, digits, dot, underscore, and hyphen.
#'
#' @param x Character vector.
#'
#' @return Character vector of the same length with unsafe characters replaced.
#'
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

#' Internal: safe maximum (NA-robust)
#'
#' Returns \code{max(x, na.rm = TRUE)}. If all values are \code{NA}, returns
#' \code{NA_real_}.
#'
#' @param x Numeric vector.
#'
#' @return Numeric scalar maximum, or \code{NA_real_}.
#'
#' @keywords internal
safe_max <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

#' Internal: safe mean (NA-robust)
#'
#' Returns \code{mean(x, na.rm = TRUE)}. If all values are \code{NA}, returns
#' \code{NA_real_}.
#'
#' @param x Numeric vector.
#'
#' @return Numeric scalar mean, or \code{NA_real_}.
#'
#' @keywords internal
safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}


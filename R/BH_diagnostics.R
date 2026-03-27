#' BH diagnostics at a headline alpha (threshold, discoveries, boundary support)
#'
#' Computes the Benjamini--Hochberg (BH) rejection threshold \eqn{t_\alpha}, the
#' number of discoveries \eqn{R_\alpha}, and a simple "boundary support" measure:
#' the number of p-values just above \eqn{t_\alpha}. Boundary support is intended
#' to indicate how sensitive the discovery set may be to small perturbations of
#' p-values near the cutoff.
#'
#' @param x A \code{dfdr_tbl} (or data frame) containing columns \code{id} and
#'   \code{p}.
#' @param alpha Numeric BH FDR level in \eqn{(0,1)}.
#' @param boundary Character. One of \code{"mult"} (multiplicative window) or
#'   \code{"add"} (additive window).
#' @param win Numeric. Relative width for the multiplicative boundary window
#'   (default 0.2). Used when \code{boundary = "mult"}, defining the interval
#'   \eqn{(t_\alpha, (1+\mathrm{win})t_\alpha]}.
#' @param delta Numeric. Additive width for the additive boundary window. Required
#'   when \code{boundary = "add"}, defining the interval \eqn{(t_\alpha, t_\alpha+\delta]}.
#'
#' @return
#' A list with two elements:
#' \describe{
#'   \item{headline}{A one-row \link[tibble:tibble]{tibble} containing BH summary
#'   diagnostics, including \code{t_alpha} (the BH threshold), \code{R_alpha} (the
#'   number of rejected hypotheses / discoveries), and \code{N_boundary} (the
#'   number of p-values in a right-neighborhood above the cutoff as determined by
#'   \code{boundary}, \code{win}, and \code{delta}).}
#'   \item{accepted}{A character vector of \code{id} values corresponding to
#'   discoveries (rows with \code{p <= t_alpha}). If no finite BH threshold exists,
#'   this is \code{character(0)}.}
#' }
#'
#' @examples
#' library(tibble)
#'
#' set.seed(1)
#' n <- 5000
#' x <- tibble(
#'   id = as.character(seq_len(n)),
#'   p = c(stats::runif(4500), stats::rbeta(500, 0.3, 1))
#' )
#'
#' out <- dfdr_bh_diagnostics(x, alpha = 0.01, boundary = "mult", win = 0.2)
#' out$headline
#' length(out$accepted)
#'
#' @export
dfdr_bh_diagnostics <- function(x,
                                alpha = 0.01,
                                boundary = c("mult", "add"),
                                win = 0.2,
                                delta = NULL) {
  boundary <- match.arg(boundary)
  x <- .safe_df(x)
  .check_has_cols(x, c("id", "p"))

  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    rlang::abort("`alpha` must be a single number in (0,1).")
  }
  if (boundary == "mult") {
    if (!is.numeric(win) || length(win) != 1 || !is.finite(win) || win <= 0) {
      rlang::abort("`win` must be a single positive number for boundary='mult'.")
    }
  } else {
    if (is.null(delta) || !is.numeric(delta) || length(delta) != 1 || !is.finite(delta) || delta <= 0) {
      rlang::abort("`delta` must be a single positive number for boundary='add'.")
    }
  }

  ok <- is.finite(x$p)
  p <- x$p[ok]
  ids <- x$id[ok]
  m <- length(p)

  t_alpha <- .bh_threshold(p, alpha)

  if (!is.finite(t_alpha)) {
    headline <- tibble::tibble(
      alpha = alpha,
      m = m,
      t_alpha = NA_real_,
      R_alpha = 0L,
      N_boundary = NA_integer_,
      boundary_type = boundary,
      win = if (boundary == "mult") win else NA_real_,
      delta = if (boundary == "add") delta else NA_real_
    )
    return(list(headline = headline, accepted = character()))
  }

  accepted_ids <- ids[p <= t_alpha]
  R_alpha <- length(accepted_ids)

  if (boundary == "mult") {
    upper <- min(1, (1 + win) * t_alpha)
  } else {
    upper <- min(1, t_alpha + delta)
  }
  N_boundary <- sum(p > t_alpha & p <= upper)

  headline <- tibble::tibble(
    alpha = alpha,
    m = m,
    t_alpha = t_alpha,
    R_alpha = as.integer(R_alpha),
    N_boundary = as.integer(N_boundary),
    boundary_type = boundary,
    win = if (boundary == "mult") win else NA_real_,
    delta = if (boundary == "add") delta else NA_real_
  )

  list(headline = headline, accepted = accepted_ids)
}

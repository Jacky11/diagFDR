#' BH diagnostics at a headline alpha (threshold, discoveries, boundary support)
#'
#' Computes the BH rejection threshold t_alpha, the number of discoveries R_alpha,
#' and "boundary support" near t_alpha (counts of p-values just above the cutoff).
#'
#' @param x A dfdr_tbl containing columns \code{id} and \code{p}.
#' @param alpha BH FDR level in (0,1).
#' @param boundary One of \code{"mult"} or \code{"add"}.
#' @param win Relative width for multiplicative boundary window (default 0.2).
#' @param delta Additive width for additive boundary window (required if boundary="add").
#'
#' @return A list with elements \code{headline} (tibble) and \code{accepted} (character vector of IDs).
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

#' Headline stability diagnostics at a given alpha
#' @param x dfdr_tbl
#' @param alpha Numeric FDR threshold between 0 and 1.
#' @param k_fixed integer for +/-K sensitivity (default 10)
#' @param k_sqrt_mult numeric multiplier for adaptive +/- ceil(mult*sqrt(D))
#' @export
dfdr_headline <- function(x, alpha = 0.01, k_fixed = 10, k_sqrt_mult = 2) {
  validate_dfdr_tbl(x)

  Ta <- sum(!x$is_decoy & x$q <= alpha)
  Da <- sum( x$is_decoy & x$q <= alpha)

  k2 <- ceiling(k_sqrt_mult * sqrt(max(Da, 0)))

  tibble::tibble(
    alpha = alpha,
    T_alpha = as.integer(Ta),
    D_alpha = as.integer(Da),
    FDR_hat = if (Ta > 0) Da / Ta else NA_real_,
    CV_hat  = if (Da > 0) 1 / sqrt(Da) else Inf,
    FDR_minus1 = if (Ta > 0) max(Da - 1, 0) / Ta else NA_real_,
    FDR_plus1  = if (Ta > 0) (Da + 1) / Ta else NA_real_,
    FDR_minusK = if (Ta > 0) max(Da - k_fixed, 0) / Ta else NA_real_,
    FDR_plusK  = if (Ta > 0) (Da + k_fixed) / Ta else NA_real_,
    k2sqrtD = as.integer(k2),
    FDR_minus2sqrtD = if (Ta > 0) max(Da - k2, 0) / Ta else NA_real_,
    FDR_plus2sqrtD  = if (Ta > 0) (Da + k2) / Ta else NA_real_
  )
}

#' Stability curve across alpha values
#' @param x An \code{dfdr_tbl}.
#' @param alphas Numeric vector of thresholds.
#' @param k_fixed Integer. Fixed +/-K sensitivity interval.
#' @param k_sqrt_mult Numeric. Multiplier for adaptive interval.
#' @return A tibble of stability metrics for each alpha.
#' @export
dfdr_curve_stability <- function(x, alphas,
                                  k_fixed = 10, k_sqrt_mult = 2) {
  validate_dfdr_tbl(x)
  purrr::map_dfr(alphas, ~dfdr_headline(x, alpha = .x, k_fixed = k_fixed, k_sqrt_mult = k_sqrt_mult))
}

#' Local tail support near the cutoff
#'
#' Counts decoys in the interval \eqn{(\alpha, \alpha + \rho \alpha]} where \eqn{\rho}
#' is \code{win_rel}. This approximates how well supported the decision boundary is.
#'
#' @param x An \code{dfdr_tbl}.
#' @param alphas Numeric vector of thresholds.
#' @param win_rel Relative window width \eqn{\rho} (default 0.2).
#' @param truncation How to handle cases where \eqn{\alpha + \rho\alpha} exceeds the
#'   available q-value export range. One of \code{"warn_drop"}, \code{"drop"}, \code{"cap"}.
#'
#' @return A tibble with \code{D_alpha_win} and related columns for each \code{alpha}.
#' @export
dfdr_local_tail <- function(x, alphas, win_rel = 0.2,
                             truncation = c("warn_drop", "drop", "cap")) {
  validate_dfdr_tbl(x)
  truncation <- match.arg(truncation)

  meta <- attr(x, "meta") %||% list()
  qmax <- meta$q_max_export %||% NA_real_
  if (!is.finite(qmax)) qmax <- max(x$q, na.rm = TRUE)

  purrr::map_dfr(alphas, function(a) {
    delta <- win_rel * a
    upper <- a + delta

    if (upper <= qmax + 1e-12) {
      win <- x[x$q > a & x$q <= upper, , drop = FALSE]
      return(tibble::tibble(
        alpha = a, delta = delta,
        n_win = nrow(win),
        D_alpha_win = sum(win$is_decoy),
        decoy_frac_win = if (nrow(win) > 0) mean(win$is_decoy) else NA_real_,
        upper_used = upper,
        truncated = FALSE
      ))
    }

    if (truncation %in% c("warn_drop", "drop")) {
      if (truncation == "warn_drop") {
        rlang::warn(paste0(
          "Local tail window exceeds q_max_export (alpha=", a,
          ", upper=", signif(upper, 4),
          ", qmax=", signif(qmax, 4), "). Returning NA."
        ))
      }
      return(tibble::tibble(
        alpha = a, delta = delta,
        n_win = NA_integer_,
        D_alpha_win = NA_integer_,
        decoy_frac_win = NA_real_,
        upper_used = NA_real_,
        truncated = TRUE
      ))
    }

    # cap
    upper2 <- qmax
    win <- x[x$q > a & x$q <= upper2, , drop = FALSE]
    tibble::tibble(
      alpha = a, delta = upper2 - a,
      n_win = nrow(win),
      D_alpha_win = sum(win$is_decoy),
      decoy_frac_win = if (nrow(win) > 0) mean(win$is_decoy) else NA_real_,
      upper_used = upper2,
      truncated = TRUE
    )
  })
}

#' Threshold elasticity (list stability to cutoff perturbation)
#'
#' Computes \eqn{S_\alpha(\epsilon) = J(A(\alpha), A((1+\epsilon)\alpha))} where
#' \eqn{J} is the Jaccard index on accepted target IDs.
#'
#' @param x An \code{dfdr_tbl}.
#' @param alphas Numeric vector of thresholds.
#' @param eps Relative perturbation (e.g. 0.2 compares alpha vs 1.2*alpha).
#'
#' @return A tibble with Jaccard overlap for each \code{alpha}.
#' @export
dfdr_elasticity <- function(x, alphas, eps = 0.2) {
  validate_dfdr_tbl(x)

  purrr::map_dfr(alphas, function(a) {
    A1 <- unique(x$id[!x$is_decoy & x$q <= a])
    A2 <- unique(x$id[!x$is_decoy & x$q <= (1 + eps) * a])
    tibble::tibble(
      alpha = a, eps = eps,
      n1 = length(A1), n2 = length(A2),
      jaccard = jaccard_vec(A1, A2)
    )
  })
}

#' Inter-run Jaccard overlap matrix at a fixed threshold
#'
#' Computes a run-by-run Jaccard similarity matrix based on accepted target IDs.
#'
#' @param x An \code{dfdr_tbl} with a non-missing \code{run} column.
#' @param alpha Threshold used for acceptance.
#'
#' @return A tibble with columns \code{run1}, \code{run2}, \code{jaccard}.
#' @export
dfdr_run_jaccard <- function(x, alpha = 0.01) {
  validate_dfdr_tbl(x)
  if (all(is.na(x$run))) rlang::abort("`run` is all NA: run-level Jaccard requires a run column.")

  dd <- x |>
    dplyr::filter(!is_decoy, q <= alpha, !is.na(run)) |>
    dplyr::select(run, id) |>
    dplyr::distinct()

  runs <- sort(unique(dd$run))
  sets <- split(dd$id, dd$run)

  out <- expand.grid(run1 = runs, run2 = runs, stringsAsFactors = FALSE)
  out$jaccard <- purrr::map2_dbl(out$run1, out$run2, function(r1, r2) {
    if (r1 == r2) return(1)
    jaccard_vec(sets[[r1]] %||% character(), sets[[r2]] %||% character())
  })
  tibble::as_tibble(out)
}

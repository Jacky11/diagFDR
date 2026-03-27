#' Headline stability diagnostics at a given alpha
#'
#' Computes headline counts at the acceptance threshold corresponding to
#' \code{alpha} (e.g., numbers of accepted targets/decoys) and simple stability
#' diagnostics that quantify how sensitive the result is to small changes in the
#' decoy boundary support.
#'
#' @param x A \code{dfdr_tbl} containing at least \code{id}, \code{is_decoy}, and a
#'   q-value column (see Details).
#' @param alpha Numeric FDR threshold in \eqn{(0,1]} at which headline metrics are
#'   computed.
#' @param k_fixed Integer. Fixed perturbation size used for sensitivity analyses
#'   based on \eqn{D \pm k} around the decoy count at the threshold.
#' @param k_sqrt_mult Numeric. Multiplier for adaptive perturbations of size
#'   \eqn{\lceil \mathrm{mult}\sqrt{D} \rceil}.
#'
#' @details
#' The function assumes that \code{x} contains target/decoy labels in
#' \code{is_decoy} and q-values in column \code{q}. If your input uses another
#' name (e.g., \code{q_value}), rename it to \code{q} before calling this
#' function.
#'
#' @return
#' A one-row \link[tibble:tibble]{tibble} with headline diagnostics evaluated at
#' the specified \code{alpha}. The table includes the number of accepted targets
#' and decoys at the threshold and additional sensitivity/stability summaries
#' derived from perturbing the boundary decoy support (using \code{k_fixed} and
#' \code{k_sqrt_mult}). The returned tibble is intended for reporting and for
#' comparison across workflows.
#'
#' @examples
#' library(tibble)
#'
#' set.seed(1)
#' n <- 5000
#' x <- tibble(
#'   id = as.character(seq_len(n)),
#'   run = "run1",
#'   is_decoy = sample(c(FALSE, TRUE), n, replace = TRUE, prob = c(0.95, 0.05)),
#'   score = rnorm(n),
#'   pep = NA_real_
#' )
#'
#' # Construct simple TDC q-values from the score (higher is better)
#' x$q <- diagFDR:::tdc_qvalues(score = x$score, is_decoy = x$is_decoy, add_decoy = 1L)
#' x <- as_dfdr_tbl(x, unit = "psm", scope = "global", q_source = "toy_tdc")
#'
#' dfdr_headline(x, alpha = 0.01)
#'
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
#'
#' Evaluates \code{\link{dfdr_headline}} over a grid of FDR thresholds
#' (\code{alphas}) to obtain a stability curve that can be plotted or compared
#' across workflows.
#'
#' @param x An \code{dfdr_tbl}.
#' @param alphas Numeric vector of FDR thresholds in \eqn{(0,1]}.
#' @param k_fixed Integer. Fixed +/-K sensitivity interval used by
#'   \code{\link{dfdr_headline}}.
#' @param k_sqrt_mult Numeric. Multiplier for the adaptive +/- ceiling(mult*sqrt(D))
#'   interval used by \code{\link{dfdr_headline}}.
#'
#' @return
#' A \link[tibble:tibble]{tibble} with one row per value of \code{alphas}.
#' Columns are those returned by \code{\link{dfdr_headline}} plus the
#' corresponding \code{alpha}. The table summarizes how headline quantities
#' (e.g., accepted target/decoy counts and sensitivity diagnostics) vary with the
#' chosen FDR cutoff.
#'
#' @examples
#' library(tibble)
#'
#' set.seed(1)
#' n <- 5000
#' df <- tibble(
#'   id = as.character(seq_len(n)),
#'   run = "run1",
#'   is_decoy = sample(c(FALSE, TRUE), n, replace = TRUE, prob = c(0.95, 0.05)),
#'   score = rnorm(n),
#'   pep = NA_real_
#' )
#'
#' df$q <- diagFDR:::tdc_qvalues(df$score, df$is_decoy, add_decoy = 1L)
#' x <- as_dfdr_tbl(df, unit = "psm", scope = "global", q_source = "toy_tdc")
#'
#' dfdr_curve_stability(x, alphas = c(1e-3, 5e-3, 1e-2))
#'
#' @export
dfdr_curve_stability <- function(x, alphas, k_fixed = 10, k_sqrt_mult = 2) {
  validate_dfdr_tbl(x)
  purrr::map_dfr(alphas, ~dfdr_headline(x, alpha = .x, k_fixed = k_fixed, k_sqrt_mult = k_sqrt_mult))
}


#' Local tail support near the cutoff
#'
#' Counts the number of observations (and decoys) in a right-neighborhood of each
#' threshold \eqn{\alpha}, defined as \eqn{(\alpha,\, \alpha + \rho\alpha]} where
#' \eqn{\rho} is \code{win_rel}. This approximates how well supported the decision
#' boundary is by decoys in the immediate tail.
#'
#' @param x An \code{dfdr_tbl}.
#' @param alphas Numeric vector of FDR thresholds in \eqn{(0,1]}.
#' @param win_rel Numeric. Relative window width \eqn{\rho} (default 0.2).
#' @param truncation Character. How to handle cases where
#'   \eqn{\alpha + \rho\alpha} exceeds the available q-value export range
#'   (i.e. \code{q_max_export}). One of \code{"warn_drop"}, \code{"drop"},
#'   \code{"cap"}.
#'
#' @return
#' A \link[tibble:tibble]{tibble} with one row per \code{alpha}. Columns include:
#' \describe{
#'   \item{alpha}{The threshold.}
#'   \item{delta}{The effective window width used (\eqn{\rho\alpha} or capped).}
#'   \item{n_win}{Number of entries with \code{q} in the window.}
#'   \item{D_alpha_win}{Number of decoys in the window.}
#'   \item{decoy_frac_win}{Decoy fraction in the window (\code{D_alpha_win/n_win}).}
#'   \item{upper_used}{Upper endpoint actually used for the window.}
#'   \item{truncated}{Logical indicating whether the requested window was truncated
#'     (or dropped) due to export limits.}
#' }
#' This table is used to assess boundary support: very small \code{n_win} or
#' \code{D_alpha_win} indicates an unstable/poorly-supported cutoff.
#'
#' @examples
#' library(tibble)
#'
#' set.seed(1)
#' n <- 5000
#' df <- tibble(
#'   id = as.character(seq_len(n)),
#'   run = "run1",
#'   is_decoy = sample(c(FALSE, TRUE), n, replace = TRUE, prob = c(0.95, 0.05)),
#'   score = rnorm(n),
#'   pep = NA_real_
#' )
#' df$q <- diagFDR:::tdc_qvalues(df$score, df$is_decoy, add_decoy = 1L)
#' x <- as_dfdr_tbl(df, unit = "psm", scope = "global", q_source = "toy_tdc")
#'
#' dfdr_local_tail(x, alphas = c(1e-3, 5e-3, 1e-2), win_rel = 0.2, truncation = "cap")
#'
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
#' Computes the stability of the accepted target set under a multiplicative
#' perturbation of the threshold:
#' \deqn{S_\alpha(\epsilon) = J(A(\alpha), A((1+\epsilon)\alpha)),}
#' where \eqn{J} is the Jaccard index and \eqn{A(\alpha)} is the set of accepted
#' target IDs at threshold \eqn{\alpha}.
#'
#' @param x An \code{dfdr_tbl}.
#' @param alphas Numeric vector of FDR thresholds in \eqn{(0,1]}.
#' @param eps Numeric scalar \eqn{\epsilon \ge 0}. Relative perturbation (e.g.
#'   \code{0.2} compares \eqn{\alpha} vs \eqn{1.2\alpha}).
#'
#' @return
#' A \link[tibble:tibble]{tibble} with one row per \code{alpha}. Columns include:
#' \describe{
#'   \item{alpha}{The baseline threshold.}
#'   \item{eps}{The relative perturbation.}
#'   \item{n1}{Number of accepted targets at \code{alpha}.}
#'   \item{n2}{Number of accepted targets at \code{(1+eps)*alpha}.}
#'   \item{jaccard}{Jaccard overlap between the two accepted target sets.}
#' }
#'
#' @examples
#' library(tibble)
#'
#' set.seed(1)
#' n <- 5000
#' df <- tibble(
#'   id = as.character(seq_len(n)),
#'   run = "run1",
#'   is_decoy = sample(c(FALSE, TRUE), n, replace = TRUE, prob = c(0.95, 0.05)),
#'   score = rnorm(n),
#'   pep = NA_real_
#' )
#'
#' df$q <- diagFDR:::tdc_qvalues(df$score, df$is_decoy, add_decoy = 1L)
#' x <- as_dfdr_tbl(df, unit = "psm", scope = "global", q_source = "toy_tdc")
#'
#' dfdr_elasticity(x, alphas = c(1e-3, 5e-3, 1e-2), eps = 0.2)
#'
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
#' Computes pairwise Jaccard similarity between runs based on accepted target IDs
#' at a fixed threshold \code{alpha}. This can reveal run-to-run disagreement or
#' heterogeneity in identifications.
#'
#' @param x An \code{dfdr_tbl} with a non-missing \code{run} column.
#' @param alpha Numeric FDR threshold in \eqn{(0,1]} used to define accepted
#'   targets (\code{q <= alpha}).
#'
#' @return
#' A \link[tibble:tibble]{tibble} with one row per run pair and columns:
#' \describe{
#'   \item{run1}{Name/identifier of the first run.}
#'   \item{run2}{Name/identifier of the second run.}
#'   \item{jaccard}{Jaccard similarity between accepted target ID sets.}
#' }
#'
#' @examples
#' library(tibble)
#'
#' set.seed(1)
#' n <- 8000
#' df <- tibble(
#'   id = as.character(seq_len(n)),
#'   run = sample(c("runA", "runB", "runC"), n, replace = TRUE),
#'   is_decoy = sample(c(FALSE, TRUE), n, replace = TRUE, prob = c(0.95, 0.05)),
#'   score = rnorm(n),
#'   pep = NA_real_
#' )
#' df$q <- diagFDR:::tdc_qvalues(df$score, df$is_decoy, add_decoy = 1L)
#' x <- as_dfdr_tbl(df, unit = "psm", scope = "global", q_source = "toy_tdc")
#'
#' dfdr_run_jaccard(x, alpha = 0.01)
#'
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

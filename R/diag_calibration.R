
#' Posterior Error Probability (PEP) reliability and internal PEP calibration error (IPE)
#'
#' Bins identifications by predicted PEP and compares mean predicted PEP to the
#' observed decoy fraction within each bin (internal target-decoy consistency check).
#'
#' @param x An \code{dfdr_tbl} with a non-missing \code{pep} column.
#' @param binwidth Bin width for PEP binning (default 0.05).
#' @param n_min Minimum bin size to include bins in the IPE summary.
#' @param pep_max Maximum mean PEP to include in the IPE summary (default 0.5).
#'
#' @return A list with elements \code{bins} (tibble) and \code{IPE} (numeric).
#' @export
dfdr_pep_reliability <- function(x, binwidth = 0.05, n_min = 200, pep_max = 0.5) {
  validate_dfdr_tbl(x)
  if (!"pep" %in% names(x) || all(is.na(x$pep))) rlang::abort("No PEP available (pep is all NA).")

  brks <- seq(0, 1, by = binwidth)
  if (utils::tail(brks, 1) < 1) brks <- c(brks, 1)

  bins <- x |>
    dplyr::filter(is.finite(.data$pep), .data$pep >= 0, .data$pep <= 1) |>
    dplyr::mutate(bin = cut(.data$pep, breaks = brks, include.lowest = TRUE)) |>
    dplyr::group_by(.data$bin) |>
    dplyr::summarise(
      pep_mean = mean(.data$pep),
      decoy_rate = mean(.data$is_decoy),
      n = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::filter(is.finite(.data$pep_mean))

  bins_use <- bins |>
    dplyr::filter(.data$n >= n_min, .data$pep_mean <= pep_max)

  IPE <- if (nrow(bins_use) > 0) {
    sum((bins_use$n / sum(bins_use$n)) * abs(bins_use$pep_mean - bins_use$decoy_rate))
  } else {
    NA_real_
  }

  list(
    bins = bins,
    IPE = IPE,
    params = list(binwidth = binwidth, n_min = n_min, pep_max = pep_max)
  )
}

#' Expected number of false targets among accepted identifications
#'
#' Computes \eqn{\sum_{i \in A(\alpha)} PEP_i} among accepted targets for each threshold.
#'
#' @param x An \code{dfdr_tbl} with non-missing \code{pep}.
#' @param alphas Numeric vector of thresholds.
#'
#' @return A tibble with \code{sum_PEP} and \code{mean_PEP} for each \code{alpha}.
#' @export
dfdr_sumpep <- function(x, alphas) {
  validate_dfdr_tbl(x)
  if (all(is.na(x$pep))) rlang::abort("No PEP available (pep is all NA).")

  purrr::map_dfr(alphas, function(a) {
    acc <- x |>
      dplyr::filter(!is_decoy, q <= a, is.finite(pep), pep >= 0, pep <= 1)
    tibble::tibble(
      alpha = a,
      n_targets = nrow(acc),
      sum_PEP = sum(acc$pep, na.rm = TRUE),
      mean_PEP = mean(acc$pep, na.rm = TRUE)
    )
  })
}

#' Equal-chance plausibility by q-value bands
#'
#' Computes decoy fractions in q-value bands and performs a pooled binomial test
#' of whether the decoy fraction is near 0.5 in a mismatch-dominated region.
#'
#' @param x An \code{dfdr_tbl}.
#' @param breaks Numeric vector of q-value cutpoints used to form bands.
#' @param low_conf Length-2 numeric vector giving the pooled test interval
#'   (e.g. \code{c(0.2, 0.5)}).
#' @param min_N Minimum pooled sample size required for the pooled test to be
#'   considered stable.
#'
#' @return A list with elements \code{bands} (tibble) and \code{pooled} (tibble).
#' @export
dfdr_equal_chance_qbands <- function(
    x,
    breaks = c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5),
    low_conf = c(0.2, 0.5),
    min_N = 2000
) {
  validate_dfdr_tbl(x)
  meta <- attr(x, "meta") %||% list()
  qmax <- meta$q_max_export %||% NA_real_
  if (!is.finite(qmax)) qmax <- max(x$q, na.rm = TRUE)

  breaks_use <- breaks[breaks <= qmax + 1e-12]
  breaks_use <- unique(sort(breaks_use))

  # --- NEW: graceful handling when not applicable ---
  if (length(breaks_use) < 2) {
    bands <- tibble::tibble(
      qbin = factor(),
      n = integer(),
      n_decoy = integer(),
      decoy_frac = numeric(),
      q_mean = numeric()
    )
    pooled <- tibble::tibble(
      qmax_export = qmax,
      low_lo = low_conf[1],
      low_hi = low_conf[2],
      N_test = 0,
      N_D_test = 0,
      pi_D_hat = NA_real_,
      effect_abs = NA_real_,
      ci95_lo = NA_real_,
      ci95_hi = NA_real_,
      p_value_binom = NA_real_,
      pass_minN = FALSE,
      note = "Not applicable: qmax below smallest nonzero break; provide smaller breaks or export higher q."
    )
    return(list(
      bands = bands,
      pooled = pooled,
      params = list(breaks = breaks, low_conf = low_conf, min_N = min_N)
    ))
  }
  # --- end NEW ---

  if (max(breaks_use) < qmax) breaks_use <- c(breaks_use, qmax)

  dd <- dplyr::mutate(
    dplyr::filter(x, is.finite(.data$q), .data$q >= 0, .data$q <= qmax),
    qbin = cut(.data$q, breaks = breaks_use, include.lowest = TRUE, right = TRUE)
  )

  bands <- dplyr::arrange(
    dplyr::summarise(
      dplyr::group_by(dd, .data$qbin),
      n = dplyr::n(),
      n_decoy = sum(.data$is_decoy),
      decoy_frac = .data$n_decoy / .data$n,
      q_mean = mean(.data$q),
      .groups = "drop"
    ),
    .data$q_mean
  )

  lo <- low_conf[1]; hi <- low_conf[2]
  test_tbl <- dplyr::filter(bands, is.finite(.data$q_mean), .data$q_mean >= lo, .data$q_mean <= hi)
  N_test <- sum(test_tbl$n)
  N_D_test <- sum(test_tbl$n_decoy)
  pi_hat <- ifelse(N_test > 0, N_D_test / N_test, NA_real_)
  ci <- wilson_ci(N_D_test, N_test, conf = 0.95)
  p_binom <- if (N_test > 0) stats::binom.test(N_D_test, N_test, p = 0.5)$p.value else NA_real_

  pooled <- tibble::tibble(
    qmax_export = qmax, low_lo = lo, low_hi = hi,
    N_test = N_test, N_D_test = N_D_test,
    pi_D_hat = pi_hat,
    effect_abs = ifelse(is.finite(pi_hat), abs(pi_hat - 0.5), NA_real_),
    ci95_lo = ci[1], ci95_hi = ci[2],
    p_value_binom = p_binom,
    pass_minN = N_test >= min_N
  )

  list(bands = bands, pooled = pooled, params = list(breaks = breaks, low_conf = low_conf, min_N = min_N))
}

#' Decoy PEP sanity checks
#'
#' Reports how many decoys receive surprisingly small PEP values.
#'
#' @param x A dfdr_tbl with a non-missing pep column.
#' @param thresholds Numeric vector of PEP thresholds to summarize.
#' @return Tibble with counts and fractions for decoys (and targets optionally).
#' @export
dfdr_pep_decoy_sanity <- function(x, thresholds = c(0.01, 0.05, 0.1, 0.2)) {
  validate_dfdr_tbl(x)
  if (!"pep" %in% names(x)) rlang::abort("`x` must contain a `pep` column.")

  pep <- x$pep
  ok <- is.finite(pep) & pep >= 0 & pep <= 1
  pep <- pep[ok]
  is_decoy <- x$is_decoy[ok]

  n_decoy <- sum(is_decoy)
  n_target <- sum(!is_decoy)

  purrr::map_dfr(thresholds, function(th) {
    decoy_le <- sum(is_decoy & pep <= th)
    target_le <- sum((!is_decoy) & pep <= th)

    tibble::tibble(
      threshold = th,
      n_decoy = n_decoy,
      n_target = n_target,
      decoy_le = decoy_le,
      target_le = target_le,
      frac_decoy_le = decoy_le / pmax(n_decoy, 1),
      frac_target_le = target_le / pmax(n_target, 1)
    )
  })
}


#' Target-focused PEP reliability using a TDC-style error proxy
#'
#' Bins targets by predicted PEP and estimates the fraction of false targets in each bin
#' using decoys as a proxy (D/T).
#'
#' @param x A dfdr_tbl with pep and is_decoy.
#' @param breaks Numeric vector of PEP bin edges.
#' @param add_decoy Additive correction to D (default 0). Use 1 for conservative small-sample correction.
#' @return Tibble with per-bin summaries.
#' @export
dfdr_pep_reliability_tdc <- function(x,
                                     breaks = seq(0, 0.5, by = 0.05),
                                     add_decoy = 0L) {
  validate_dfdr_tbl(x)
  if (!"pep" %in% names(x)) rlang::abort("`x` must contain a `pep` column.")

  breaks <- sort(unique(breaks))
  if (length(breaks) < 2) rlang::abort("`breaks` must have at least two values.")

  df <- x |>
    dplyr::filter(is.finite(.data$pep), .data$pep >= min(breaks), .data$pep <= max(breaks)) |>
    dplyr::mutate(
      bin_i = findInterval(.data$pep, vec = breaks, rightmost.closed = TRUE),
      bin_i = dplyr::if_else(.data$bin_i < 1L | .data$bin_i >= length(breaks), NA_integer_, .data$bin_i)
    ) |>
    dplyr::filter(!is.na(.data$bin_i))

  mids <- (breaks[-length(breaks)] + breaks[-1]) / 2

  out <- df |>
    dplyr::group_by(.data$bin_i) |>
    dplyr::summarise(
      pep_lo = breaks[unique(.data$bin_i)],
      pep_hi = breaks[unique(.data$bin_i) + 1L],
      pep_bin_mid = mids[unique(.data$bin_i)],
      n_target = sum(!.data$is_decoy),
      n_decoy = sum(.data$is_decoy),
      pep_mean_targets = mean(.data$pep[!.data$is_decoy], na.rm = TRUE),
      pep_mean_all = mean(.data$pep, na.rm = TRUE),
      err_hat_target = pmin(1, (.data$n_decoy + add_decoy) / pmax(.data$n_target, 1)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      bin = sprintf("(%.3g, %.3g]", .data$pep_lo, .data$pep_hi)
    ) |>
    dplyr::select(.data$bin, .data$n_target, .data$n_decoy,
                  .data$pep_mean_targets, .data$pep_mean_all,
                  .data$err_hat_target, .data$pep_bin_mid)

  out
}

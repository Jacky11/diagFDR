#' Flag diagnostic values based on heuristic thresholds
#'
#' Adds simple, ASCII-only flags and a coarse overall status to a headline
#' diagnostics table (typically produced by \code{\link{dfdr_headline}} or the
#' \code{headline} table inside \code{\link{dfdr_run_all}}). The heuristics are
#' intended for quick triage in reports; they do not replace manual review.
#'
#' Rendering as icons (if desired) can be handled downstream (e.g. in an HTML
#' template); this function returns plain text labels for portability.
#'
#' @param headline_tbl A data.frame/tibble of headline diagnostics. It may contain
#'   one or multiple rows. If some expected columns are missing, they are created
#'   and filled with \code{NA}.
#'
#' @return
#' A \link[tibble:tibble]{tibble} (or data frame) with the same rows as
#' \code{headline_tbl} and additional columns:
#' \describe{
#'   \item{flag_Dalpha}{Flag based on \code{D_alpha} (decoy support).}
#'   \item{flag_CV}{Flag based on \code{CV_hat} (stability/variability).}
#'   \item{flag_Dwin}{Flag based on \code{D_alpha_win} (local tail support).}
#'   \item{flag_IPE}{Flag based on \code{IPE} (internal PEP calibration error).}
#'   \item{flag_FDR}{Flag comparing \code{FDR_hat} to the nominal \code{alpha}.}
#'   \item{flag_equalchance}{Flag based on \code{effect_abs} (equal-chance deviation).}
#'   \item{status}{Overall triage status: \code{"VALID"}, \code{"CAUTION"}, or
#'   \code{"REVIEW_NEEDED"}.}
#'   \item{interpretation}{A short human-readable interpretation string.}
#' }
#'
#' @examples
#' library(tibble)
#'
#' headline_tbl <- tibble(
#'   list = c("A", "B"),
#'   alpha = 0.01,
#'   D_alpha = c(5, 80),
#'   CV_hat = c(0.35, 0.10),
#'   D_alpha_win = c(2, 30),
#'   IPE = c(0.12, 0.02),
#'   FDR_hat = c(0.02, 0.009),
#'   effect_abs = c(0.18, 0.05)
#' )
#'
#' flag_headline(headline_tbl)
#'
#' @export
flag_headline <- function(headline_tbl) {
  if (!is.data.frame(headline_tbl)) {
    rlang::abort("`headline_tbl` must be a data.frame/tibble.")
  }

  # Ensure required columns exist (some diagnostics are optional)
  ensure_col <- function(df, nm, default = NA_real_) {
    if (!nm %in% names(df)) df[[nm]] <- default
    df
  }

  headline_tbl <- headline_tbl |>
    ensure_col("D_alpha", NA_real_) |>
    ensure_col("CV_hat", NA_real_) |>
    ensure_col("D_alpha_win", NA_real_) |>
    ensure_col("IPE", NA_real_) |>
    ensure_col("FDR_hat", NA_real_) |>
    ensure_col("alpha", NA_real_) |>
    ensure_col("effect_abs", NA_real_)

  headline_tbl |>
    dplyr::mutate(
      flag_Dalpha = dplyr::case_when(
        is.na(D_alpha) ~ "NA",
        D_alpha < 10 ~ "PROBLEM_granular",
        D_alpha < 50 ~ "CAUTION",
        TRUE ~ "OK"
      ),
      flag_CV = dplyr::case_when(
        is.na(CV_hat) ~ "NA",
        CV_hat > 0.30 ~ "PROBLEM_granular",
        CV_hat > 0.15 ~ "CAUTION",
        TRUE ~ "OK"
      ),
      flag_Dwin = dplyr::case_when(
        is.na(D_alpha_win) ~ "NA",
        D_alpha_win < 5 ~ "PROBLEM_very_sparse",
        D_alpha_win < 20 ~ "CAUTION_sparse",
        TRUE ~ "OK"
      ),
      flag_IPE = dplyr::case_when(
        is.na(IPE) ~ "NA",
        IPE > 0.10 ~ "PROBLEM_serious",
        IPE > 0.05 ~ "CAUTION",
        TRUE ~ "OK"
      ),
      flag_FDR = dplyr::case_when(
        is.na(FDR_hat) ~ "NA",
        is.na(alpha) ~ "NA",
        FDR_hat > alpha * 1.5 ~ "PROBLEM_inflated",
        FDR_hat > alpha * 1.2 ~ "CAUTION_slightly_high",
        TRUE ~ "OK"
      ),
      flag_equalchance = dplyr::case_when(
        is.na(effect_abs) ~ "NA",
        effect_abs > 0.15 ~ "PROBLEM_imbalance",
        effect_abs > 0.10 ~ "CAUTION",
        TRUE ~ "OK"
      ),

      status = dplyr::case_when(
        grepl(
          "PROBLEM",
          paste(flag_Dalpha, flag_CV, flag_Dwin, flag_IPE, flag_FDR, flag_equalchance)
        ) ~ "REVIEW_NEEDED",
        grepl(
          "CAUTION",
          paste(flag_Dalpha, flag_CV, flag_Dwin, flag_IPE, flag_FDR, flag_equalchance)
        ) ~ "CAUTION",
        TRUE ~ "VALID"
      ),

      interpretation = dplyr::case_when(
        flag_FDR == "PROBLEM_inflated" ~ "Empirical FDR exceeds nominal threshold - likely scope misuse",
        flag_Dalpha == "PROBLEM_granular" & flag_CV == "PROBLEM_granular" ~ "Sparse decoy support - threshold potentially unstable",
        flag_IPE == "PROBLEM_serious" ~ "PEPs seriously miscalibrated",
        flag_equalchance == "PROBLEM_imbalance" ~ "Equal-chance assumption potentially violated - check if selection effects were used (two-pass learning, etc.)",
        status == "VALID" ~ "No major issues detected",
        TRUE ~ "Minor issues - review flagged metrics"
      )
    )
}

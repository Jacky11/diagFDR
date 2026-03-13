#' Flag diagnostic values based on heuristic thresholds
#'
#' Produces ASCII-only flags for portability. Rendering as icons is handled in the HTML template.
#'
#' @param headline_tbl Output from dfdr_headline (possibly stacked)
#' @return Same table with added columns: status, flags, interpretation
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

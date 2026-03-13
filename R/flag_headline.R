#' Flag diagnostic values based on heuristic thresholds
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
        is.na(D_alpha) ~ "n/a",
        D_alpha < 10 ~ "⚠️ granular",
        D_alpha < 50 ~ "⚠ caution",
        TRUE ~ "✓"
      ),
      flag_CV = dplyr::case_when(
        is.na(CV_hat) ~ "n/a",
        CV_hat > 0.30 ~ "⚠️ granular",
        CV_hat > 0.15 ~ "⚠ caution",
        TRUE ~ "✓"
      ),
      flag_Dwin = dplyr::case_when(
        is.na(D_alpha_win) ~ "n/a",
        D_alpha_win < 5 ~ "⚠️ very sparse",
        D_alpha_win < 20 ~ "⚠ sparse",
        TRUE ~ "✓"
      ),
      flag_IPE = dplyr::case_when(
        is.na(IPE) ~ "n/a",
        IPE > 0.10 ~ "⚠️ serious",
        IPE > 0.05 ~ "⚠ caution",
        TRUE ~ "✓"
      ),
      flag_FDR = dplyr::case_when(
        is.na(FDR_hat) ~ "n/a",
        is.na(alpha) ~ "n/a",
        FDR_hat > alpha * 1.5 ~ "⚠️ inflated",
        FDR_hat > alpha * 1.2 ~ "⚠ slightly high",
        TRUE ~ "✓"
      ),
      flag_equalchance = dplyr::case_when(
        is.na(effect_abs) ~ "n/a",
        effect_abs > 0.15 ~ "⚠️ imbalance",
        effect_abs > 0.10 ~ "⚠ caution",
        TRUE ~ "✓"
      ),

      status = dplyr::case_when(
        grepl(
          "⚠️",
          paste(flag_Dalpha, flag_CV, flag_Dwin, flag_IPE, flag_FDR, flag_equalchance)
        ) ~ "⚠️ Review needed",
        grepl(
          "⚠",
          paste(flag_Dalpha, flag_CV, flag_Dwin, flag_IPE, flag_FDR, flag_equalchance)
        ) ~ "⚠ Caution",
        TRUE ~ "✓ Valid"
      ),

      interpretation = dplyr::case_when(
        flag_FDR == "⚠️ inflated" ~ "Empirical FDR exceeds nominal threshold - likely scope misuse",
        flag_Dalpha == "⚠️ granular" & flag_CV == "⚠️ granular" ~ "Sparse decoy support - threshold potentially unstable",
        flag_IPE == "⚠️ serious" ~ "PEPs seriously miscalibrated",
        flag_equalchance == "⚠️ imbalance" ~ "Equal-chance assumption potentially violated - check if selection effects has been used (two-pass learning, etc.)",
        status == "✓ Valid" ~ "No major issues detected",
        TRUE ~ "Minor issues - review flagged metrics"
      )
    )
}

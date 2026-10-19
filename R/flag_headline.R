#' Flag headline diagnostics as review alerts
#'
#' Adds ASCII flags and a coarse triage status to a headline diagnostics table
#' (typically the \code{headline} table of \code{\link{dfdr_run_all}}). The flags
#' are \emph{review alerts} on internal (target-decoy) diagnostics. They are not a
#' pass/fail verdict: a list that raises no alert has passed internal checks
#' only, and decoy validity itself requires external evidence (e.g. entrapment).
#'
#' Thresholds follow the diagFDR paper (guidelines, not cut-offs). Their
#' evidential basis differs (derived check, simulation-informed alert, heuristic
#' alert; see \code{\link{dfdr_alert_basis}}):
#' \itemize{
#'   \item decoy count: alert if \eqn{D_\alpha < 10}, caution if \eqn{< 50};
#'   \item decoy-support index \eqn{D_\alpha^{-1/2}}: alert if \eqn{> 0.30},
#'     caution if \eqn{> 0.15};
#'   \item local boundary support \eqn{D_{\alpha,win}}: alert if \eqn{< 5},
#'     caution if \eqn{< 20};
#'   \item decoy-reference PEP concordance, IPE (target-decoy comparator):
#'     alert if \eqn{> 0.05} (heuristic);
#'   \item recomputed decoy ratio: alert if \eqn{D_\alpha/T_\alpha > 1.5\alpha},
#'     caution if \eqn{> 1.2\alpha};
#'   \item equal chance under paired target-decoy competition (decoy fraction
#'     in the low-confidence q band): alert if
#'     \eqn{|\hat\pi_D - 0.5| > 0.15}, caution if \eqn{> 0.10}; reported as not
#'     evaluable (\code{"NE_two_pass"}) when the list comes from a two-pass
#'     export, whose low-confidence region is conditioned by the first pass. For a
#'     DIA-NN list whose pass status was not given (\code{two_pass = NA}), an
#'     imbalance is reported as \code{"CAUTION_imbalance_pass_unknown"}.
#' }
#'
#' @param headline_tbl A data.frame/tibble of headline diagnostics. Missing
#'   expected columns are created and filled with \code{NA}. Optional columns
#'   \code{support_index}, \code{eq_evaluable} and \code{eq_reason} are used when
#'   present.
#'
#' @return
#' A \link[tibble:tibble]{tibble} with the rows of \code{headline_tbl} and
#' additional columns \code{flag_Dalpha}, \code{flag_support} (decoy-support
#' index), \code{flag_Dwin}, \code{flag_IPE}, \code{flag_FDR} (recomputed
#' decoy ratio), \code{flag_equalchance}, \code{status} (\code{"NO_ALERT"},
#' \code{"CAUTION"} or \code{"ALERT"}) and \code{interpretation}. Flag values are
#' \code{"OK"}, \code{"CAUTION*"}, \code{"ALERT_*"}, \code{"NE_*"} (not
#' evaluable) or \code{"NA"} (not computed).
#'
#' @examples
#' library(tibble)
#'
#' headline_tbl <- tibble(
#'   list = c("A", "B"),
#'   alpha = 0.01,
#'   D_alpha = c(5, 80),
#'   support_index = 1 / sqrt(c(5, 80)),
#'   D_alpha_win = c(2, 30),
#'   IPE = c(0.12, 0.02),
#'   FDR_hat = c(0.02, 0.009),
#'   effect_abs = c(0.18, 0.05),
#'   eq_evaluable = c(TRUE, FALSE),
#'   eq_reason = c(NA, "two_pass")
#' )
#'
#' flag_headline(headline_tbl)
#'
#' @export
flag_headline <- function(headline_tbl) {
  if (!is.data.frame(headline_tbl)) {
    rlang::abort("`headline_tbl` must be a data.frame/tibble.")
  }

  ensure_col <- function(df, nm, default = NA_real_) {
    if (!nm %in% names(df)) df[[nm]] <- default
    df
  }

  headline_tbl <- headline_tbl |>
    ensure_col("D_alpha", NA_real_) |>
    
    ensure_col("support_index", NA_real_) |>
    ensure_col("D_alpha_win", NA_real_) |>
    ensure_col("IPE", NA_real_) |>
    ensure_col("FDR_hat", NA_real_) |>
    ensure_col("alpha", NA_real_) |>
    ensure_col("effect_abs", NA_real_) |>
    ensure_col("eq_evaluable", NA) |>
    ensure_col("eq_reason", NA_character_)

    # support index from D_alpha when the column is missing
  si <- ifelse(is.na(headline_tbl$support_index) & is.finite(headline_tbl$D_alpha) & headline_tbl$D_alpha > 0,
               1 / sqrt(headline_tbl$D_alpha), headline_tbl$support_index)
  headline_tbl$support_index <- si


  out <- headline_tbl |>
    dplyr::mutate(
      flag_Dalpha = dplyr::case_when(
        is.na(D_alpha) ~ "NA",
        D_alpha < 10 ~ "ALERT_granular",
        D_alpha < 50 ~ "CAUTION",
        TRUE ~ "OK"
      ),
      flag_support = dplyr::case_when(
        is.na(support_index) ~ "NA",
        support_index > 0.30 ~ "ALERT_granular",
        support_index > 0.15 ~ "CAUTION",
        TRUE ~ "OK"
      ),
      
      flag_Dwin = dplyr::case_when(
        is.na(D_alpha_win) ~ "NA",
        D_alpha_win < 5 ~ "ALERT_sparse",
        D_alpha_win < 20 ~ "CAUTION_sparse",
        TRUE ~ "OK"
      ),
      flag_IPE = dplyr::case_when(
        is.na(IPE) ~ "NA",
        IPE > 0.05 ~ "ALERT_above_0.05",
        TRUE ~ "OK"
      ),
      flag_FDR = dplyr::case_when(
        is.na(FDR_hat) ~ "NA",
        is.na(alpha) ~ "NA",
        FDR_hat > alpha * 1.5 ~ "ALERT_ratio_above_alpha",
        FDR_hat > alpha * 1.2 ~ "CAUTION_ratio_above_alpha",
        TRUE ~ "OK"
      ),
      flag_equalchance = dplyr::case_when(
        !is.na(eq_reason) & eq_reason == "two_pass" ~ "NE_two_pass",
        !is.na(eq_evaluable) & !eq_evaluable ~ "NE_insufficient_N",
        is.na(effect_abs) ~ "NA",
                effect_abs > 0.15 & !is.na(eq_reason) & eq_reason == "pass_unknown" ~
          "CAUTION_imbalance_pass_unknown",
        effect_abs > 0.15 ~ "ALERT_imbalance",
        effect_abs > 0.10 ~ "CAUTION",
        TRUE ~ "OK"
      )
    )

  all_flags <- paste(out$flag_Dalpha, out$flag_support, out$flag_Dwin,
                     out$flag_IPE, out$flag_FDR, out$flag_equalchance)

  out |>
    dplyr::mutate(
      status = dplyr::case_when(
        grepl("ALERT", all_flags) ~ "ALERT",
        grepl("CAUTION", all_flags) ~ "CAUTION",
        TRUE ~ "NO_ALERT"
      ),
      interpretation = dplyr::case_when(
        flag_FDR == "ALERT_ratio_above_alpha" ~
          "Recomputed decoy ratio D/T exceeds alpha: check which q-value column was thresholded and how rows were aggregated (scope)",
        flag_Dalpha == "ALERT_granular" | flag_support == "ALERT_granular" ~
          "Few accepted decoys: the realised FDP is poorly constrained at this threshold",
        flag_IPE == "ALERT_above_0.05" ~
          "PEPs disagree with the decoy-derived reference (check aggregation of PEPs across runs)",
        flag_equalchance == "ALERT_imbalance" ~
          "Decoy fraction at low confidence departs from 0.5: equal chance under paired competition questionable (or selection by an earlier pass)",
                flag_equalchance == "CAUTION_imbalance_pass_unknown" ~
          "Decoy fraction at low confidence departs from 0.5; if this export is two-pass (DIA-NN --reanalyse), set two_pass = TRUE: the check is then not evaluable",
        status == "NO_ALERT" ~ "No internal alert (internal checks only; decoy validity needs external evidence)",
        TRUE ~ "Minor deviations - review flagged metrics"
      )
    )
}

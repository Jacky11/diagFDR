#' Build a headline summary table (one row per list)
#'
#' @param diag Output of dfdr_run_all().
#' @return Tibble with one row per list, suitable for exporting as summary_headline.csv.
#' @export
dfdr_summary_headline <- function(diag) {
  stopifnot(is.list(diag))
  head_tbl <- diag$tables$headline
  if (!is.data.frame(head_tbl) || nrow(head_tbl) == 0) return(tibble::tibble())

  alpha_main <- (diag$params %||% list())$alpha_main %||% NA_real_

  out <- head_tbl

  # ---- local boundary support at alpha_main: D_alpha_win, n_win, decoy_frac_win, truncated
  lt <- diag$tables$local_tail
  if (is.data.frame(lt) && nrow(lt) > 0 && is.finite(alpha_main)) {
    lt_at <- dplyr::filter(lt, alpha == alpha_main) |>
      dplyr::select(dplyr::any_of(c("list","D_alpha_win","n_win","decoy_frac_win","delta","upper_used","truncated")))
    out <- dplyr::left_join(out, lt_at, by = "list")
  }

  # ---- equal-chance pooled (per list already in your outputs)
  ec <- diag$tables$equal_chance_pooled
  if (is.data.frame(ec) && nrow(ec) > 0) {
    # Your table does not currently include list; dfdr_run_all mutates list=.y, so it should.
    # Keep the key columns you showed:
    keep_ec <- intersect(names(ec), c("list","low_lo","low_hi","N_test","pi_D_hat","effect_abs",
                                      "ci95_lo","ci95_hi","p_value_binom","pass_minN","qmax_export"))
    if ("list" %in% keep_ec) out <- dplyr::left_join(out, ec[, keep_ec, drop = FALSE], by = "list")
  }

  # ---- IPE (if present)
  pep_ipe <- diag$tables$pep_IPE
  if (is.data.frame(pep_ipe) && nrow(pep_ipe) > 0) {
    out <- dplyr::left_join(out, pep_ipe, by = "list")
  }

  # Keep a stable column order (adjust as you like)
  preferred <- c(
    "list","alpha",
    "T_alpha","D_alpha","FDR_hat","CV_hat",
    "D_alpha_win","n_win","decoy_frac_win","truncated",
    "FDR_minus1","FDR_plus1","FDR_minus2sqrtD","FDR_plus2sqrtD",
    "pi_D_hat","effect_abs","ci95_lo","ci95_hi","p_value_binom","pass_minN",
    "IPE"
  )
  out |>
    dplyr::select(dplyr::any_of(preferred), dplyr::everything()) |>
    dplyr::arrange(list)
}


#' Build a headline summary table (one row per list)
#'
#' Convenience helper to assemble a compact, export-ready summary table from the
#' output of \code{\link{dfdr_run_all}}. The function merges headline metrics with
#' selected columns from other diagnostics (e.g. local tail support, equal-chance
#' pooled estimate, and PEP reliability headline if available).
#'
#' @param diag A list as returned by \code{\link{dfdr_run_all}}.
#'
#' @return
#' A \link[tibble:tibble]{tibble} with one row per element of the input list used
#' in \code{\link{dfdr_run_all}} (i.e. one row per \code{list} label). The table
#' is suitable for exporting to a CSV (e.g. \code{summary_headline.csv}).
#'
#' The output typically contains (when available) headline target/decoy counts at
#' \code{alpha_main}, local boundary-support columns (e.g. \code{D_alpha_win},
#' \code{n_win}), equal-chance pooled columns (e.g. \code{pi_D_hat},
#' \code{effect_abs}), and PEP reliability headline quantities (e.g. \code{IPE}).
#' If \code{diag$tables$headline} is missing or empty, an empty tibble is returned.
#'
#' @examples
#' library(tibble)
#'
#' set.seed(1)
#' n <- 4000
#' df <- tibble(
#'   id = as.character(seq_len(n)),
#'   run = sample(c("run1", "run2"), n, replace = TRUE),
#'   is_decoy = sample(c(FALSE, TRUE), n, replace = TRUE, prob = c(0.95, 0.05)),
#'   score = rnorm(n)
#' )
#' df$q <- diagFDR:::tdc_qvalues(df$score, df$is_decoy, add_decoy = 1L)
#' df$pep <- NA_real_
#' x <- as_dfdr_tbl(df, unit = "psm", scope = "global", q_source = "toy_tdc")
#'
#' diag <- dfdr_run_all(
#'   xs = list(toy = x),
#'   alpha_main = 0.01,
#'   compute_pseudo_pvalues = FALSE
#' )
#'
#' dfdr_summary_headline(diag)
#'
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


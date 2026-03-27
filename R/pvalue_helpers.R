#' Internal: assert required columns exist
#'
#' Checks that a data frame contains the required column names.
#'
#' @param x A data.frame / tibble.
#' @param cols Character vector of required column names.
#'
#' @return
#' Returns \code{invisible(TRUE)} if all columns are present; otherwise throws an error.
#'
#' @keywords internal
.check_has_cols <- function(x, cols) {
  miss <- setdiff(cols, colnames(x))
  if (length(miss)) rlang::abort(paste0("Missing required column(s): ", paste(miss, collapse = ", ")))
}

#' Internal: validate data-frame input
#'
#' Asserts that \code{x} is a data.frame (including tibble). Used to provide a
#' consistent error message for user-facing functions that accept data frames.
#'
#' @param x Object to check.
#'
#' @return
#' Returns \code{x} unchanged if it is a data.frame; otherwise throws an error.
#'
#' @keywords internal
.safe_df <- function(x) {
  if (!inherits(x, "data.frame")) rlang::abort("`x` must be a data.frame / tibble (preferably dfdr_tbl).")
  x
}

#' Internal: extract stratification columns
#'
#' Returns the subset of \code{x} corresponding to columns named in \code{stratify}.
#' Used by stratified diagnostics.
#'
#' @param x A data.frame / tibble.
#' @param stratify Optional character vector of column names used for stratification.
#'
#' @return
#' If \code{stratify} is \code{NULL} or empty, returns \code{NULL}. Otherwise
#' returns a data.frame with the selected stratification columns (in the order
#' given by \code{stratify}). Throws an error if any requested columns are missing.
#'
#' @keywords internal
.get_strata_df <- function(x, stratify = NULL) {
  if (is.null(stratify) || length(stratify) == 0) return(NULL)
  if (!all(stratify %in% names(x))) {
    rlang::abort(paste0("Unknown `stratify` column(s): ", paste(setdiff(stratify, names(x)), collapse = ", ")))
  }
  x[stratify]
}

#' Internal: build stratum identifiers
#'
#' Builds a human-readable stratum identifier per row from a data.frame of
#' stratification columns. Each identifier is of the form
#' \code{"col1=... | col2=..."}.
#'
#' @param strata_df A data.frame of stratification columns (typically produced by
#'   \code{.get_strata_df}). If \code{NULL}, a single stratum \code{"all"} is used.
#'
#' @return
#' A character vector of stratum identifiers (length \code{nrow(strata_df)} if not
#' \code{NULL}, otherwise \code{"all"}).
#'
#' @keywords internal
.make_stratum_id <- function(strata_df) {
  if (is.null(strata_df)) return(rep("all", 1))
  apply(strata_df, 1, function(row) {
    paste(paste0(names(strata_df), "=", ifelse(is.na(row), "NA", as.character(row))), collapse = " | ")
  })
}

#' Internal: Jaccard similarity of two sets
#'
#' Computes the Jaccard index \eqn{|A \cap B| / |A \cup B|} on unique values.
#' Returns 1 when both sets are empty.
#'
#' @param a Vector representing set A.
#' @param b Vector representing set B.
#'
#' @return Numeric scalar in \eqn{[0,1]} giving the Jaccard similarity.
#'
#' @keywords internal
.jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(1)
  length(intersect(a, b)) / length(union(a, b))
}

#' Internal: Benjamini--Hochberg rejection threshold
#'
#' Computes the BH cutoff \eqn{t_\alpha} for a vector of p-values \code{p} at level
#' \code{alpha}. Returns \code{NA} if no hypotheses are rejected (or no finite p-values).
#'
#' @param p Numeric vector of p-values (only finite values are used).
#' @param alpha Numeric scalar in \eqn{(0,1)}.
#'
#' @return Numeric scalar BH threshold \eqn{t_\alpha}, or \code{NA_real_} if no
#'   rejection threshold exists.
#'
#' @keywords internal
.bh_threshold <- function(p, alpha) {
  ok <- is.finite(p)
  p <- p[ok]
  m <- length(p)
  if (m == 0) return(NA_real_)
  o <- order(p)
  p_sorted <- p[o]
  crit <- (seq_len(m) / m) * alpha
  idx <- which(p_sorted <= crit)
  if (length(idx) == 0) return(NA_real_)
  p_sorted[max(idx)]
}

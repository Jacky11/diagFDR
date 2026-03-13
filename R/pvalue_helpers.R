#' @keywords internal
.check_has_cols <- function(x, cols) {
  miss <- setdiff(cols, colnames(x))
  if (length(miss)) rlang::abort(paste0("Missing required column(s): ", paste(miss, collapse = ", ")))
}

#' @keywords internal
.safe_df <- function(x) {
  if (!inherits(x, "data.frame")) rlang::abort("`x` must be a data.frame / tibble (preferably dfdr_tbl).")
  x
}

#' @keywords internal
.get_strata_df <- function(x, stratify = NULL) {
  if (is.null(stratify) || length(stratify) == 0) return(NULL)
  if (!all(stratify %in% names(x))) {
    rlang::abort(paste0("Unknown `stratify` column(s): ", paste(setdiff(stratify, names(x)), collapse = ", ")))
  }
  x[stratify]
}

#' @keywords internal
.make_stratum_id <- function(strata_df) {
  if (is.null(strata_df)) return(rep("all", 1))
  apply(strata_df, 1, function(row) {
    paste(paste0(names(strata_df), "=", ifelse(is.na(row), "NA", as.character(row))), collapse = " | ")
  })
}

#' @keywords internal
.jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(1)
  length(intersect(a, b)) / length(union(a, b))
}

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


#' Scope disagreement between two lists
#'
#' Computes Jaccard overlap of accepted targets between two \code{dfdr_tbl} objects
#' across thresholds. Both inputs must use the same \code{id} semantics.
#'
#' @param x1 First \code{dfdr_tbl}.
#' @param x2 Second \code{dfdr_tbl}.
#' @param alphas Numeric vector of thresholds.
#' @param label1 Label for \code{x1} in output.
#' @param label2 Label for \code{x2} in output.
#' @param normalize_ids Logical; if TRUE (default), extracts base precursor ID
#'   by stripping run prefixes (everything before "||"). Set to FALSE if IDs
#'   are already in the same format.
#'
#' @return A tibble with Jaccard overlap for each \code{alpha}.
#' @export
dfdr_scope_disagreement <- function(x1, x2, alphas,
                                    label1 = "A", label2 = "B",
                                    normalize_ids = TRUE) {
  validate_dfdr_tbl(x1)
  validate_dfdr_tbl(x2)

  # Check if ID formats differ
  has_pipe_x1 <- any(grepl("\\|\\|", x1$id[1:min(100, nrow(x1))]))
  has_pipe_x2 <- any(grepl("\\|\\|", x2$id[1:min(100, nrow(x2))]))

  if (has_pipe_x1 != has_pipe_x2 && !normalize_ids) {
    rlang::warn(
      "ID formats appear to differ. Consider setting `normalize_ids = TRUE`."
    )
  }

  # Helper to extract base precursor ID
  extract_base_id <- function(ids) {
    if (!normalize_ids) return(ids)
    ifelse(grepl("\\|\\|", ids),
           sub("^.*\\|\\|", "", ids),
           ids)
  }

  purrr::map_dfr(alphas, function(a) {
    A <- unique(extract_base_id(x1$id[!x1$is_decoy & x1$q <= a]))
    B <- unique(extract_base_id(x2$id[!x2$is_decoy & x2$q <= a]))
    tibble::tibble(
      alpha = a,
      label1 = label1,
      label2 = label2,
      n1 = length(A),
      n2 = length(B),
      jaccard = jaccard_vec(A, B)
    )
  })
}

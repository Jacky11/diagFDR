#' Scope disagreement between two lists
#'
#' Computes the Jaccard overlap of accepted target IDs between two
#' \code{dfdr_tbl} objects across a set of thresholds. This is useful for comparing
#' results obtained under different FDR scopes (e.g. run-wise vs global) or from
#' different pipelines, provided that both inputs use compatible \code{id}
#' semantics.
#'
#' @param x1 First \code{dfdr_tbl}.
#' @param x2 Second \code{dfdr_tbl}.
#' @param alphas Numeric vector of thresholds in \eqn{(0,1]}.
#' @param label1 Character label for \code{x1} in the output (default \code{"A"}).
#' @param label2 Character label for \code{x2} in the output (default \code{"B"}).
#' @param normalize_ids Logical. If \code{TRUE} (default), extracts a "base" ID by
#'   stripping any run prefix of the form \code{"<run>||<id>"} (everything up to
#'   and including \code{"||"}). Set to \code{FALSE} if IDs are already comparable
#'   between \code{x1} and \code{x2}.
#'
#' @return
#' A \link[tibble:tibble]{tibble} with one row per \code{alpha} and columns:
#' \describe{
#'   \item{alpha}{The threshold.}
#'   \item{label1,label2}{The labels identifying the two inputs.}
#'   \item{n1,n2}{Numbers of accepted target IDs in \code{x1} and \code{x2} at \code{alpha}.}
#'   \item{jaccard}{Jaccard similarity between the two accepted target ID sets.}
#' }
#'
#' @examples
#' library(tibble)
#'
#' set.seed(1)
#' n <- 4000
#'
#' # Two "lists" with slightly different q-values for the same base IDs
#' base_ids <- as.character(seq_len(n))
#'
#' df1 <- tibble(
#'   id = paste0("run1||", base_ids),
#'   run = "run1",
#'   is_decoy = sample(c(FALSE, TRUE), n, replace = TRUE, prob = c(0.95, 0.05)),
#'   score = rnorm(n),
#'   q = stats::runif(n, 0, 0.05),
#'   pep = NA_real_
#' )
#'
#' df2 <- tibble(
#'   id = base_ids,              # no run prefix here
#'   run = "all",
#'   is_decoy = df1$is_decoy,    # same decoy labels for simplicity
#'   score = df1$score + rnorm(n, sd = 0.2),
#'   q = pmin(1, df1$q * 1.2),   # slightly perturbed q-values
#'   pep = NA_real_
#' )
#'
#' x1 <- as_dfdr_tbl(df1, unit = "psm", scope = "runwise", q_source = "toy")
#' x2 <- as_dfdr_tbl(df2, unit = "psm", scope = "global",  q_source = "toy")
#'
#' dfdr_scope_disagreement(x1, x2, alphas = c(0.005, 0.01, 0.02), normalize_ids = TRUE)
#'
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

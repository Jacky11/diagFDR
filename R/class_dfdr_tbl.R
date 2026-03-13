#' Validate a dfdr_tbl
#' @keywords internal
validate_dfdr_tbl <- function(x) {
  stopifnot(is.data.frame(x))

  req <- c("id", "is_decoy", "q", "pep", "run", "score")
  miss <- setdiff(req, names(x))
  if (length(miss)) {
    rlang::abort(paste0("Missing required columns: ", paste(miss, collapse = ", ")))
  }

  if (!is.character(x$id)) rlang::abort("`id` must be character.")

  if (!is.logical(x$is_decoy)) rlang::abort("`is_decoy` must be logical (TRUE/FALSE).")
  if (anyNA(x$is_decoy)) rlang::abort("`is_decoy` contains NA; please fix upstream.")

  if (!is.numeric(x$q)) rlang::abort("`q` must be numeric.")
  if (any(!is.finite(x$q))) rlang::abort("`q` contains non-finite values.")
  if (any(x$q < 0 | x$q > 1)) rlang::abort("`q` must be in [0,1].")

  if (!is.numeric(x$pep)) rlang::abort("`pep` must be numeric (use NA if unavailable).")
  if (any(is.finite(x$pep) & (x$pep < 0 | x$pep > 1))) {
    rlang::abort("`pep` must be in [0,1] when present.")
  }

  if (!is.character(x$run)) rlang::abort("`run` must be character (use NA if not applicable).")

  if (!is.numeric(x$score)) rlang::abort("`score` must be numeric (use NA if unavailable).")

  # Optional p-values / pseudo-p-values
  if ("p" %in% names(x)) {
    if (!is.numeric(x$p)) rlang::abort("`p` must be numeric (use NA if unavailable).")
    if (any(is.finite(x$p) & (x$p < 0 | x$p > 1))) {
      rlang::abort("`p` must be in [0,1] when present.")
    }
  }

  x
}

#' Create an dfdr_tbl
#'
#' @param x A data.frame with columns id,is_decoy,q,pep,run,score (pep/run/score can be NA).
#' @param unit Character. E.g. "precursor", "psm", "peptide", "protein", "run×precursor".
#' @param scope Character. E.g. "runwise", "global", "aggregated".
#' @param q_source Character describing where q came from.
#' @param q_max_export Numeric. Optional export ceiling used to generate q (if known).
#' @param p_source Optional character describing where p came from.
#' @param provenance Named list (tool, version, params, command).
#'
#' @return An object of class \code{dfdr_tbl} (inheriting from \code{tbl_df}) containing the validated input data with \code{is_decoy} coerced to logical. Metadata (\code{unit}, \code{scope}, \code{q_source}, \code{q_max_export}, \code{provenance}) are stored as attributes and accessible via \code{attr(x, "meta")}.
#' @export
as_dfdr_tbl <- function(x,
                        unit = NA_character_,
                        scope = NA_character_,
                        q_source = NA_character_,
                        q_max_export = NA_real_,
                        p_source = NA_character_,
                        provenance = list()) {
  x <- tibble::as_tibble(x)

  if (!"id" %in% names(x)) rlang::abort("Missing required column: id")
  if (!"is_decoy" %in% names(x)) rlang::abort("Missing required column: is_decoy")

  # Coerce is_decoy if needed
  if (is.numeric(x$is_decoy) || is.integer(x$is_decoy)) {
    u <- sort(unique(stats::na.omit(x$is_decoy)))
    if (!all(u %in% c(0, 1))) {
      rlang::abort("Numeric `is_decoy` must contain only 0/1 (no NAs).")
    }
    x$is_decoy <- x$is_decoy == 1
  }
  if (is.character(x$is_decoy)) {
    if (all(stats::na.omit(x$is_decoy) %in% c("0","1"))) {
      x$is_decoy <- x$is_decoy == "1"
    } else {
      x$is_decoy <- as.logical(x$is_decoy)
    }
  }

  x <- validate_dfdr_tbl(x)

  attr(x, "meta") <- list(
    unit = unit,
    scope = scope,
    q_source = q_source,
    q_max_export = q_max_export,
    p_source = p_source,
    provenance = provenance
  )
  class(x) <- c("dfdr_tbl", class(x))
  x
}

#' @export
print.dfdr_tbl <- function(x, ...) {
  meta <- attr(x, "meta") %||% list()
  cat("<dfdr_tbl>\n")
  cat("  unit:    ", meta$unit %||% NA_character_, "\n", sep = "")
  cat("  scope:   ", meta$scope %||% NA_character_, "\n", sep = "")
  cat("  q_from:  ", meta$q_source %||% NA_character_, "\n", sep = "")
  if (!is.null(meta$p_source) && !is.na(meta$p_source)) {
    cat("  p_from:  ", meta$p_source, "\n", sep = "")
  }
  cat("  n:       ", nrow(x), "\n", sep = "")
  NextMethod()
}

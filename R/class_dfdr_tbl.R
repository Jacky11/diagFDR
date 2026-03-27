#' Validate a \code{dfdr_tbl}
#'
#' Performs structural checks required by \code{diagFDR} diagnostics. This
#' function is intended for internal use and is called by constructors and
#' downstream diagnostics.
#'
#' The required columns are \code{id}, \code{is_decoy}, \code{q}, \code{pep},
#' \code{run}, and \code{score}. If a \code{p} column is present, it is also
#' validated.
#'
#' @param x A data frame (typically a tibble) to validate.
#'
#' @return
#' The input \code{x} invisibly (unchanged) if validation succeeds. Otherwise an
#' error is raised describing the first detected problem.
#'
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

#' Create a \code{dfdr_tbl}
#'
#' Coerces input to a tibble, standardizes the \code{is_decoy} column to logical,
#' validates required columns and types, attaches metadata, and returns an S3
#' object of class \code{dfdr_tbl} for downstream diagnostics.
#'
#' @param x A data.frame with columns \code{id}, \code{is_decoy}, \code{q},
#'   \code{pep}, \code{run}, \code{score} (where \code{pep}/\code{run}/\code{score}
#'   may be \code{NA} if unavailable). Optionally may contain \code{p}.
#' @param unit Character. Statistical unit, e.g. \code{"psm"}, \code{"precursor"},
#'   \code{"peptide"}, \code{"protein"}.
#' @param scope Character. Scope of FDR control, e.g. \code{"runwise"},
#'   \code{"global"}, \code{"aggregated"}.
#' @param q_source Character. Description of where q-values came from (e.g. a
#'   column name or tool).
#' @param q_max_export Numeric. Optional export ceiling used to generate q (if known).
#' @param p_source Optional character describing where \code{p} came from.
#' @param provenance Named list carrying provenance information (tool, version,
#'   parameters, command, etc.).
#'
#' @return
#' An S3 object of class \code{dfdr_tbl} that inherits from
#' \link[tibble:tibble]{tbl_df} (and \code{data.frame}). The returned object:
#' \itemize{
#'   \item contains the validated input data with \code{is_decoy} coerced to
#'   logical (\code{TRUE}/\code{FALSE});
#'   \item carries metadata in \code{attr(x, "meta")} with entries \code{unit},
#'   \code{scope}, \code{q_source}, \code{q_max_export}, \code{p_source}, and
#'   \code{provenance}.
#' }
#'
#' @examples
#' library(tibble)
#'
#' df <- tibble(
#'   id = as.character(1:6),
#'   run = c("run1","run1","run1","run2","run2","run2"),
#'   is_decoy = c(FALSE, FALSE, TRUE, FALSE, TRUE, FALSE),
#'   score = c(10, 9, 8, 7, 6, 5),
#'   q = c(0.01, 0.02, 0.02, 0.05, 0.06, 0.07),
#'   pep = c(0.01, 0.02, NA, 0.05, NA, 0.07)
#' )
#'
#' x <- as_dfdr_tbl(
#'   df,
#'   unit = "psm",
#'   scope = "global",
#'   q_source = "toy_q",
#'   provenance = list(tool = "toy")
#' )
#'
#' x
#' attr(x, "meta")
#'
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

#' Print a \code{dfdr_tbl}
#'
#' Prints a compact header with attached metadata (unit, scope, q-source, etc.)
#' and then falls back to the default tibble/data.frame print method.
#'
#' @param x A \code{dfdr_tbl}.
#' @param ... Passed to the next print method.
#'
#' @return
#' Returns \code{x} invisibly (called for its printing side effect).
#'
#' @examples
#' library(tibble)
#'
#' df <- tibble(
#'   id = as.character(1:3),
#'   run = "run1",
#'   is_decoy = c(FALSE, TRUE, FALSE),
#'   score = c(10, 9, 8),
#'   q = c(0.01, 0.02, 0.03),
#'   pep = c(0.01, NA, 0.03)
#' )
#' x <- as_dfdr_tbl(df, unit = "psm", scope = "global", q_source = "toy_q")
#' x
#'
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

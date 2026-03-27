#' Detect DIA-NN run column name
#' @keywords internal
diann_detect_run_col <- function(df, candidates = c("Run", "File.Name", "File.Name.Index")) {
  nm <- names(df)
  hit <- candidates[candidates %in% nm]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

#' Validate that required DIA-NN columns exist
#' @keywords internal
diann_require_cols <- function(df, cols) {
  miss <- setdiff(cols, names(df))
  if (length(miss)) rlang::abort(paste0("Missing required DIA-NN columns: ", paste(miss, collapse = ", ")))
  invisible(TRUE)
}

#' Read a DIA-NN \code{report.parquet}
#'
#' Reads a DIA-NN \code{report.parquet} file using the \pkg{arrow} package
#' (suggested dependency) and returns it as a tibble.
#'
#' @param path Character scalar. Path to a DIA-NN \code{report.parquet} file.
#'
#' @return
#' A \link[tibble:tibble]{tibble} containing the columns present in the DIA-NN
#' parquet report (column names depend on DIA-NN export settings).
#'
#' @examples
#' if (requireNamespace("arrow", quietly = TRUE)) {
#'   # Create a tiny parquet file and read it back
#'   tmp <- tempfile(fileext = ".parquet")
#'   df <- data.frame(Precursor.Id = c("P1", "P2"), Q.Value = c(0.01, 0.02))
#'   arrow::write_parquet(df, tmp)
#'   out <- read_diann_parquet(tmp)
#'   out
#' }
#'
#' @export
read_diann_parquet <- function(path) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    rlang::abort("Reading DIA-NN parquet requires the 'arrow' package (install.packages('arrow')).")
  }
  if (!file.exists(path)) rlang::abort(paste0("File not found: ", path))
  arrow::read_parquet(path) |> tibble::as_tibble()
}

#' DIA-NN -> global precursor universe (deduplicated by \code{Precursor.Id})
#'
#' Constructs a precursor-level universe by aggregating a DIA-NN report table to
#' one row per \code{Precursor.Id}. For each precursor, the minimum q-value is
#' used (via \code{safe_min}) and \code{is_decoy} is set to \code{TRUE} if any row
#' in the group is a decoy.
#'
#' @param rep A DIA-NN report tibble (e.g. returned by \code{\link{read_diann_parquet}}).
#' @param q_col Character. Column name for q-values (default \code{"Q.Value"}).
#' @param pep_col Optional character. PEP column name (default \code{"PEP"} if present).
#' @param score_col Optional character. Score column name to retain; if missing, \code{score} is \code{NA}.
#' @param q_max_export Optional numeric export ceiling used in DIA-NN export (e.g. 0.5).
#' @param unit Character. Unit metadata stored in the returned object.
#' @param scope Character. Scope metadata stored in the returned object.
#' @param q_source Character. Label stored in metadata describing the q-value source.
#'
#' @return
#' A \code{dfdr_tbl} (tibble subclass) with one row per precursor and required
#' columns \code{id}, \code{is_decoy}, \code{q}, \code{pep}, \code{run}, \code{score}.
#' Metadata are stored in \code{attr(x, "meta")}.
#'
#' @examples
#' library(tibble)
#'
#' rep <- tibble(
#'   Precursor.Id = c("P1", "P1", "P2"),
#'   Decoy = c(0L, 0L, 1L),
#'   Q.Value = c(0.01, 0.02, 0.03),
#'   PEP = c(0.02, 0.01, 0.9)
#' )
#'
#' x <- diann_global_precursor(rep, q_col = "Q.Value", q_max_export = 0.5)
#' x
#'
#' @export
diann_global_precursor <- function(rep,
                                   q_col = "Q.Value",
                                   pep_col = NULL,
                                   score_col = NULL,
                                   q_max_export = NA_real_,
                                   unit = "precursor",
                                   scope = "global",
                                   q_source = q_col) {
  diann_require_cols(rep, c("Precursor.Id", "Decoy", q_col))

  if (is.null(pep_col)) pep_col <- if ("PEP" %in% names(rep)) "PEP" else NULL

  # score is not consistently named across DIA-NN exports; keep optional
  has_score <- !is.null(score_col) && score_col %in% names(rep)
  has_pep <- !is.null(pep_col) && pep_col %in% names(rep)

  dd <- rep |>
    dplyr::transmute(
      id = as.character(.data[["Precursor.Id"]]),
      is_decoy = as.integer(.data[["Decoy"]]) == 1L,
      q = suppressWarnings(as.numeric(.data[[q_col]])),
      pep = if (has_pep) suppressWarnings(as.numeric(.data[[pep_col]])) else NA_real_,
      run = NA_character_,
      score = if (has_score) suppressWarnings(as.numeric(.data[[score_col]])) else NA_real_
    ) |>
    dplyr::filter(!is.na(id), is.finite(q), q >= 0, q <= 1) |>
    dplyr::group_by(id) |>
    dplyr::summarise(
      is_decoy = any(is_decoy),
      q = safe_min(q),
      pep = safe_min(pep),
      run = NA_character_,
      score = safe_min(score),
      .groups = "drop"
    )

  as_dfdr_tbl(
    dd,
    unit = unit,
    scope = scope,
    q_source = q_source,
    q_max_export = q_max_export,
    provenance = list(tool = "DIA-NN")
  )
}

#' DIA-NN -> run-by-precursor universe (deduplicated within run)
#'
#' Constructs a run-specific universe by aggregating a DIA-NN report table to one
#' row per \code{(run, Precursor.Id)}. The run column is detected automatically
#' unless provided. For each group, the minimum q-value is used and \code{is_decoy}
#' is \code{TRUE} if any row is a decoy.
#'
#' @param rep A DIA-NN report tibble.
#' @param q_col Character. q-value column name (default \code{"Q.Value"}).
#' @param run_col Optional character. Name of the run column. If \code{NULL}, the
#'   function attempts to detect one of \code{"Run"}, \code{"File.Name"}, or
#'   \code{"File.Name.Index"}.
#' @param pep_col Optional character. PEP column name (default \code{"PEP"} if present).
#' @param score_col Optional character. Score column name (if present).
#' @param q_max_export Optional numeric export ceiling (e.g. 0.5), stored as metadata.
#' @param id_mode Character. Either \code{"runxid"} (default) to set
#'   \code{id = paste(run, Precursor.Id, sep="||")} or \code{"id"} to set
#'   \code{id = Precursor.Id} and keep \code{run} separate.
#' @param unit Character. Unit metadata stored in the returned object.
#' @param scope Character. Scope metadata stored in the returned object.
#' @param q_source Character. Label stored in metadata describing the q-value source.
#'
#' @return
#' A \code{dfdr_tbl} with one row per run-by-precursor unit. The returned table
#' includes \code{id}, \code{run}, \code{is_decoy}, \code{q}, \code{pep}, and
#' \code{score}. Metadata are stored in \code{attr(x, "meta")}.
#'
#' @examples
#' library(tibble)
#'
#' rep <- tibble(
#'   Run = c("r1", "r1", "r2", "r2"),
#'   Precursor.Id = c("P1", "P1", "P1", "P2"),
#'   Decoy = c(0L, 0L, 1L, 0L),
#'   Q.Value = c(0.01, 0.02, 0.03, 0.02),
#'   PEP = c(0.02, 0.01, 0.9, 0.05)
#' )
#'
#' x <- diann_runxprecursor(rep, q_col = "Q.Value", run_col = "Run", id_mode = "runxid")
#' x
#'
#' @export
diann_runxprecursor <- function(rep,
                                q_col = "Q.Value",
                                run_col = NULL,
                                pep_col = NULL,
                                score_col = NULL,
                                q_max_export = NA_real_,
                                id_mode = c("id", "runxid"),
                                unit = "runxprecursor",
                                scope = "runwise",
                                q_source = q_col) {
  diann_require_cols(rep, c("Precursor.Id", "Decoy", q_col))

  id_mode <- match.arg(id_mode)

  if (is.null(run_col)) run_col <- diann_detect_run_col(rep)
  if (is.na(run_col)) rlang::abort("No run column found (expected one of Run / File.Name / File.Name.Index).")

  if (is.null(pep_col)) pep_col <- if ("PEP" %in% names(rep)) "PEP" else NULL
  has_pep <- !is.null(pep_col) && pep_col %in% names(rep)
  has_score <- !is.null(score_col) && score_col %in% names(rep)

  dd <- rep |>
    dplyr::transmute(
      run = as.character(.data[[run_col]]),
      precursor = as.character(.data[["Precursor.Id"]]),
      is_decoy = as.integer(.data[["Decoy"]]) == 1L,
      q = suppressWarnings(as.numeric(.data[[q_col]])),
      pep = if (has_pep) suppressWarnings(as.numeric(.data[[pep_col]])) else NA_real_,
      score = if (has_score) suppressWarnings(as.numeric(.data[[score_col]])) else NA_real_
    ) |>
    dplyr::filter(!is.na(run), !is.na(precursor), is.finite(q), q >= 0, q <= 1) |>
    dplyr::group_by(run, precursor) |>
    dplyr::summarise(
      is_decoy = any(is_decoy),
      q = safe_min(q),
      pep = safe_min(pep),
      score = safe_min(score),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      id = if (id_mode == "runxid") paste(run, precursor, sep = "||") else precursor
    ) |>
    dplyr::select(id, is_decoy, q, pep, run, score)

  as_dfdr_tbl(
    dd,
    unit = unit,
    scope = scope,
    q_source = q_source,
    q_max_export = q_max_export,
    provenance = list(tool = "DIA-NN", run_col = run_col, id_mode = id_mode)
  )
}

#' DIA-NN: global precursor list built by minimum run-wise q across runs
#'
#' Constructs a precursor-level table by taking the minimum run-wise q-value
#' across runs for each \code{Precursor.Id}. This mirrors a common scope misuse
#' (aggregating run-wise q-values into a single global list) and is useful for
#' scope-disagreement diagnostics.
#'
#' @param rep A DIA-NN report table (e.g. returned by \code{\link{read_diann_parquet}}).
#' @param run_col Optional character. Name of the run column. If \code{NULL},
#'   attempts to detect it (e.g. \code{"Run"} or \code{"File.Name"}).
#' @param q_col Character. Name of the q-value column to aggregate (default \code{"Q.Value"}).
#' @param pep_col Optional character. Name of the PEP column (default \code{"PEP"} if present).
#' @param score_col Optional character. Name of a score column to carry along (if present).
#' @param q_max_export Optional numeric export ceiling used by the tool (e.g. 0.5),
#'   stored as metadata for truncation-aware diagnostics.
#' @param unit Character. Unit metadata stored in the returned object.
#' @param scope Character. Scope metadata stored in the returned object.
#' @param q_source Character. Source label stored in metadata.
#'
#' @return
#' A \code{dfdr_tbl} with one row per precursor, where \code{q} equals the minimum
#' of the run-wise q-values across runs for that precursor. The returned object is
#' intended for diagnostics (e.g. scope disagreement), not as a recommended FDR procedure.
#'
#' @examples
#' library(tibble)
#'
#' rep <- tibble(
#'   Run = c("r1", "r2", "r1", "r2"),
#'   Precursor.Id = c("P1", "P1", "P2", "P2"),
#'   Decoy = c(0L, 0L, 0L, 1L),
#'   Q.Value = c(0.02, 0.01, 0.05, 0.04),
#'   PEP = c(0.03, 0.02, 0.1, 0.9)
#' )
#'
#' x <- diann_global_minrunq(rep, run_col = "Run", q_col = "Q.Value", q_max_export = 0.5)
#' x
#'
#' @export
diann_global_minrunq <- function(rep,
                                 run_col = NULL,
                                 q_col = "Q.Value",
                                 pep_col = NULL,
                                 score_col = NULL,
                                 q_max_export = NA_real_,
                                 unit = "precursor",
                                 scope = "aggregated",
                                 q_source = paste0("min_run(", q_col, ")")) {
  diann_require_cols(rep, c("Precursor.Id", "Decoy", q_col))

  if (is.null(run_col)) run_col <- diann_detect_run_col(rep)
  if (is.na(run_col)) rlang::abort("No run column found (expected one of Run / File.Name / File.Name.Index).")

  if (is.null(pep_col)) pep_col <- if ("PEP" %in% names(rep)) "PEP" else NULL
  has_pep <- !is.null(pep_col) && pep_col %in% names(rep)
  has_score <- !is.null(score_col) && score_col %in% names(rep)

  dd <- rep |>
    dplyr::transmute(
      run = as.character(.data[[run_col]]),
      id = as.character(.data[["Precursor.Id"]]),
      is_decoy = as.integer(.data[["Decoy"]]) == 1L,
      q = suppressWarnings(as.numeric(.data[[q_col]])),
      pep = if (has_pep) suppressWarnings(as.numeric(.data[[pep_col]])) else NA_real_,
      score = if (has_score) suppressWarnings(as.numeric(.data[[score_col]])) else NA_real_
    ) |>
    dplyr::filter(!is.na(run), !is.na(id), is.finite(q), q >= 0, q <= 1) |>
    dplyr::group_by(id) |>
    dplyr::summarise(
      is_decoy = any(is_decoy),
      q = safe_min(q),     # min over runs
      pep = safe_min(pep), # min PEP over runs (not guaranteed calibrated; still useful as diagnostic)
      run = NA_character_,
      score = safe_min(score),
      .groups = "drop"
    ) |>
    dplyr::select(id, is_decoy, q, pep, run, score)

  as_dfdr_tbl(
    dd,
    unit = unit,
    scope = scope,
    q_source = q_source,
    q_max_export = q_max_export,
    provenance = list(tool = "DIA-NN", run_col = run_col)
  )
}

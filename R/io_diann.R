#' Read a DIA-NN report.parquet
#'
#' Requires the 'arrow' package (in Suggests).
#'
#' @param path Path to report.parquet
#' @return tibble with DIA-NN columns
#' @export
read_diann_parquet <- function(path) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    rlang::abort("Reading DIA-NN parquet requires the 'arrow' package (install.packages('arrow')).")
  }
  if (!file.exists(path)) rlang::abort(paste0("File not found: ", path))
  arrow::read_parquet(path) |> tibble::as_tibble()
}

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

#' DIA-NN -> global precursor universe (deduplicated by Precursor.Id)
#'
#' Uses min q per precursor (safe_min), and marks Decoy TRUE if any row in group is decoy.
#'
#' @param rep DIA-NN report tibble (from read_diann_parquet)
#' @param q_col Column name for q-values (default "Q.Value" or "Global.Q.Value")
#' @param pep_col Optional PEP column name (default "PEP" if present)
#' @param score_col Optional score column name (if you want to store it; otherwise NA)
#' @param q_max_export Optional export ceiling used in DIA-NN export (e.g. 0.5)
#' @param unit unit metadata
#' @param scope scope metadata
#' @param q_source q_source metadata
#' @return dfdr_tbl
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

#' DIA-NN -> runxprecursor universe (deduplicated within run by Precursor.Id)
#'
#' @param rep DIA-NN report tibble
#' @param q_col q-value column name (default "Q.Value")
#' @param run_col optional; auto-detect if NULL
#' @param pep_col optional PEP column
#' @param score_col optional score column
#' @param q_max_export export ceiling (e.g. 0.5)
#' @return dfdr_tbl with run column filled and id = paste(run, Precursor.Id, sep="||") OR
#'         id = Precursor.Id depending on id_mode.
#' @param id_mode "runxid" (default) sets id = run||precursor; "id" sets id=Precursor.Id and keeps run separate.
#' @param unit Unit metadata stored in the returned object.
#' @param scope Scope metadata stored in the returned object.
#' @param q_source Source label stored in metadata.
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

#' DIA-NN: global precursor list built by min run-wise q across runs
#'
#' Constructs a precursor-level table by taking the minimum run-wise q-value
#' (e.g. \code{Q.Value}) across runs for each \code{Precursor.Id}. This mirrors a
#' common scope misuse and is useful for scope-disagreement diagnostics.
#'
#' @param rep A DIA-NN report table (e.g. returned by \code{\link{read_diann_parquet}}).
#' @param run_col Optional. Name of the run column. If \code{NULL}, attempts to detect it
#'   (e.g. \code{"Run"} or \code{"File.Name"}).
#' @param q_col Name of the q-value column to aggregate (default \code{"Q.Value"}).
#' @param pep_col Optional. Name of the PEP column (default \code{"PEP"} if present).
#' @param score_col Optional. Name of a score column to carry along (if present).
#' @param q_max_export Optional numeric export ceiling used by the tool (e.g. 0.5).
#'   Stored as metadata for truncation-aware diagnostics.
#' @param unit Unit metadata stored in the returned object (default \code{"precursor"}).
#' @param scope Scope metadata stored in the returned object (default \code{"aggregated"}).
#' @param q_source Source label stored in metadata (default \code{"min_run(Q.Value)"}).
#'
#' @return An \code{dfdr_tbl}.
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

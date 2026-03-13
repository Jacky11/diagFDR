#' Read a Spectronaut TSV report
#'
#' Uses data.table::fread for efficient reading of large files.
#' Spectronaut uses comma as decimal separator in some locales;
#' adjust dec argument if needed.
#'
#' @param path Path to Spectronaut report (TSV)
#' @param dec Decimal separator (default ".")
#' @param select Optional character vector of columns to read (NULL = all)
#' @return data.table with Spectronaut columns
#' @export
read_spectronaut <- function(path, dec = ".", select = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    rlang::abort("Reading Spectronaut requires 'data.table' package (install.packages('data.table')).")
  }
  if (!file.exists(path)) rlang::abort(paste0("File not found: ", path))

  # Spectronaut exports can be huge; fread is much faster than read.delim
  dt <- data.table::fread(
    path,
    sep = "\t",
    dec = dec,
    select = select,
    showProgress = TRUE
  )

  tibble::as_tibble(dt)
}

#' Read Spectronaut efficiently with column selection
#'
#' For very large files, read only the columns needed for diagnostics.
#'
#' @param path Path to Spectronaut TSV
#' @param minimal If TRUE, read only essential columns for precursor-level diagnostics
#' @param dec Decimal separator (use "," for European locales)
#' @return tibble
#' @export
read_spectronaut_efficient <- function(path, minimal = TRUE, dec = ".") {
  if (minimal) {
    select_cols <- c(
      "R.FileName", "R.Condition",
      "EG.PrecursorId", "EG.IsDecoy",
      "EG.Qvalue", "EG.GlobalPrecursorQvalue",
      "EG.PEP", "EG.Cscore", "EG.NormalizedCscore"
    )
    read_spectronaut(path, dec = dec, select = select_cols)
  } else {
    read_spectronaut(path, dec = dec, select = NULL)
  }
}

#' Detect Spectronaut run column name
#' @keywords internal
spectronaut_detect_run_col <- function(df, candidates = c("R.FileName", "R.Condition", "File.Name")) {
  nm <- names(df)
  hit <- candidates[candidates %in% nm]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

#' Validate that required Spectronaut columns exist
#' @keywords internal
spectronaut_require_cols <- function(df, cols) {
  miss <- setdiff(cols, names(df))
  if (length(miss)) rlang::abort(paste0("Missing required Spectronaut columns: ", paste(miss, collapse = ", ")))
  invisible(TRUE)
}

#' Parse Spectronaut boolean decoy column
#' @keywords internal
spectronaut_parse_decoy <- function(x) {
  # Spectronaut uses "True"/"False" or TRUE/FALSE
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "1")
}

#' Recompute q-values from scores using target-decoy competition
#'
#' When Spectronaut exports NaN q-values for decoys, we can recompute them
#' from scores. Higher scores are better (more confident).
#'
#' @param df Data frame with columns: is_decoy, score
#' @return Data frame with added column: q_from_score
#' @keywords internal
spectronaut_recompute_q_from_score <- function(df) {
  df |>
    dplyr::arrange(dplyr::desc(.data$score)) |>
    dplyr::mutate(
      rank = dplyr::row_number(),
      cumsum_decoy = cumsum(.data$is_decoy),
      cumsum_total = dplyr::row_number(),
      q_from_score = pmin(1, .data$cumsum_decoy / .data$cumsum_total)
    ) |>
    dplyr::arrange(dplyr::desc(.data$rank)) |>
    dplyr::mutate(q_from_score = cummin(.data$q_from_score)) |>
    dplyr::arrange(.data$rank) |>
    dplyr::select(-.data$rank, -.data$cumsum_decoy, -.data$cumsum_total)
}

#' Spectronaut -> runxprecursor universe (deduplicated within run by EG.PrecursorId)
#'
#'
#' @param rep Spectronaut report tibble (from read_spectronaut or read_spectronaut_efficient)
#' @param q_col q-value column name (default "EG.Qvalue" for run-wise FDR control)
#' @param run_col optional; auto-detect if NULL (looks for R.FileName or R.Condition)
#' @param pep_col optional PEP column (default "EG.PEP")
#' @param score_col optional score column (default "EG.Cscore")
#' @param q_max_export export ceiling (default 1.0)
#' @param recompute_q If TRUE and decoys lack valid q-values, recompute from score per run.
#'   This is typically necessary as Spectronaut exports NaN q-values for decoys.
#' @param id_mode "id" (default) sets id=EG.PrecursorId and keeps run separate;
#'   "runxid" sets id = run||precursor for global uniqueness
#' @param unit Unit metadata (default "runxprecursor")
#' @param scope Scope metadata (default "runwise")
#' @param q_source Source label for metadata
#' @return dfdr_tbl with run column filled, suitable for run-wise FDR diagnostics
#'
#' @details
#' When Spectronaut exports decoys with NaN q-values (typical behavior), this function
#' recomputes q-values from scores per run using target-decoy competition. This respects
#' the statistical structure where each run×precursor is an independent hypothesis and
#' avoids pooling scores that may not be on comparable scales across runs.
#'
#'
#' @examples
#' \dontrun{
#' # Read Spectronaut report (use dec="," for European exports)
#' rep <- read_spectronaut_efficient("spectronaut_report.tsv", dec = ",")
#'
#' # Create run-wise universe for FDR diagnostics
#' univ <- spectronaut_runxprecursor(rep, q_col = "EG.Qvalue")
#'
#' # Run diagnostics
#' results <- dfdr_run_all(list(runwise = univ), derive_p = FALSE)
#' }
#'
#' @export
spectronaut_runxprecursor <- function(rep,
                                      q_col = "EG.Qvalue",
                                      run_col = NULL,
                                      pep_col = "EG.PEP",
                                      score_col = "EG.Cscore",
                                      q_max_export = 1.0,
                                      recompute_q = TRUE,
                                      id_mode = c("id", "runxid"),
                                      unit = "runxprecursor",
                                      scope = "runwise",
                                      q_source = q_col) {
  spectronaut_require_cols(rep, c("EG.PrecursorId", "EG.IsDecoy"))

  id_mode <- match.arg(id_mode)

  if (is.null(run_col)) run_col <- spectronaut_detect_run_col(rep)
  if (is.na(run_col)) rlang::abort("No run column found (expected R.FileName or R.Condition).")

  # Check if we need to recompute q-values
  needs_recompute <- FALSE
  if (recompute_q && !is.null(score_col) && score_col %in% names(rep)) {
    n_decoy <- sum(rep$EG.IsDecoy == TRUE)
    n_decoy_finite_q <- sum(rep$EG.IsDecoy == TRUE & is.finite(rep[[q_col]]))
    if (n_decoy > 0 && n_decoy_finite_q == 0) {
      needs_recompute <- TRUE
      message("Spectronaut decoys have NaN q-values; recomputing PER RUN from ", score_col)
    }
  }

  has_pep <- !is.null(pep_col) && pep_col %in% names(rep)
  has_score <- !is.null(score_col) && score_col %in% names(rep)

  dd <- rep |>
    dplyr::transmute(
      run = as.character(.data[[run_col]]),
      precursor = as.character(.data[["EG.PrecursorId"]]),
      is_decoy = spectronaut_parse_decoy(.data[["EG.IsDecoy"]]),
      q_original = if (q_col %in% names(rep)) suppressWarnings(as.numeric(.data[[q_col]])) else NA_real_,
      pep = if (has_pep) suppressWarnings(as.numeric(.data[[pep_col]])) else NA_real_,
      score = if (has_score) suppressWarnings(as.numeric(.data[[score_col]])) else NA_real_
    ) |>
    dplyr::filter(!is.na(run), !is.na(precursor))

  # Recompute q-values from scores if needed (PER RUN to respect run-level hypotheses)
  if (needs_recompute) {
    dd <- dd |>
      dplyr::filter(is.finite(score)) |>
      dplyr::group_by(run) |>
      dplyr::group_modify(~ spectronaut_recompute_q_from_score(.x)) |>
      dplyr::ungroup() |>
      dplyr::mutate(q = .data$q_from_score) |>
      dplyr::select(-.data$q_from_score, -(.data$q_original))
    q_source <- paste0(q_source, " (recomputed per run from ", score_col, ")")
  } else {
    dd <- dd |>
      dplyr::mutate(q = .data$q_original) |>
      dplyr::select(-.data$q_original) |>
      dplyr::filter(is.finite(q), q >= 0, q <= 1)
  }

  # Deduplicate within run (in case of duplicates within same run)
  dd <- dd |>
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
    provenance = list(tool = "Spectronaut", run_col = run_col, id_mode = id_mode, recomputed = needs_recompute)
  )
}


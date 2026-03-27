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

#' Read a Spectronaut TSV report
#'
#' Reads a Spectronaut report exported as tab-separated values (TSV). Uses
#' \code{data.table::fread()} for efficient reading of large files.
#'
#' Spectronaut may use a comma as decimal separator in some locales; set
#' \code{dec=","} when needed.
#'
#' @param path Character scalar. Path to a Spectronaut report TSV file.
#' @param dec Character scalar. Decimal separator passed to \code{fread()} (default \code{"."}).
#' @param select Optional character vector of column names to read; \code{NULL} reads all columns.
#'
#' @return
#' A \link[tibble:tibble]{tibble} containing the columns present in the Spectronaut
#' report (column names depend on the export configuration).
#'
#' @examples
#' if (requireNamespace("data.table", quietly = TRUE)) {
#'   # Create a tiny TSV and read it back
#'   tmp <- tempfile(fileext = ".tsv")
#'   txt <- paste0(
#'     "R.FileName\tEG.PrecursorId\tEG.IsDecoy\tEG.Qvalue\tEG.PEP\tEG.Cscore\n",
#'     "run1\tP1\tFalse\t0.01\t0.02\t120\n",
#'     "run1\tP2\tTrue\tNaN\t0.90\t10\n"
#'   )
#'   writeLines(txt, tmp)
#'   read_spectronaut(tmp)
#' }
#'
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
#' Convenience wrapper around \code{\link{read_spectronaut}}. For very large
#' Spectronaut reports, reading only a minimal set of columns can substantially
#' speed up I/O and reduce memory usage.
#'
#' @param path Character scalar. Path to a Spectronaut TSV report.
#' @param minimal Logical. If \code{TRUE} (default), reads only columns typically
#'   needed for precursor-level diagnostics; if \code{FALSE}, reads all columns.
#' @param dec Character scalar. Decimal separator (use \code{","} for some locales).
#'
#' @return
#' A \link[tibble:tibble]{tibble} containing either the selected minimal columns
#' (if \code{minimal=TRUE}) or all columns (if \code{minimal=FALSE}).
#'
#' @examples
#' if (requireNamespace("data.table", quietly = TRUE)) {
#'   tmp <- tempfile(fileext = ".tsv")
#'   txt <- paste0(
#'     "R.FileName\tR.Condition\tEG.PrecursorId\tEG.IsDecoy\tEG.Qvalue\tEG.PEP\tEG.Cscore\n",
#'     "run1\tcondA\tP1\tFalse\t0.01\t0.02\t120\n",
#'     "run1\tcondA\tP2\tTrue\tNaN\t0.90\t10\n"
#'   )
#'   writeLines(txt, tmp)
#'   read_spectronaut_efficient(tmp, minimal = TRUE)
#' }
#'
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

#' Spectronaut -> run-by-precursor universe (deduplicated within run by \code{EG.PrecursorId})
#'
#' Constructs a run-wise precursor universe from a Spectronaut report. The output
#' is suitable for run-wise FDR diagnostics because each row corresponds to a
#' \code{run} \eqn{\times} \code{precursor} hypothesis.
#'
#' If Spectronaut exports \code{NaN} q-values for decoys (common), the function can
#' optionally recompute q-values from scores per run using target-decoy competition.
#'
#' @param rep Spectronaut report tibble (e.g. returned by \code{\link{read_spectronaut}}
#'   or \code{\link{read_spectronaut_efficient}}).
#' @param q_col Character. q-value column name (default \code{"EG.Qvalue"} for run-wise control).
#' @param run_col Optional character. Run column name; if \code{NULL}, attempts to
#'   detect one of \code{"R.FileName"} or \code{"R.Condition"}.
#' @param pep_col Optional character. PEP column name (default \code{"EG.PEP"}).
#' @param score_col Optional character. Score column name (default \code{"EG.Cscore"}).
#' @param q_max_export Numeric. Export ceiling for q-values (default 1.0), stored in metadata.
#' @param recompute_q Logical. If \code{TRUE} (default) and decoys lack finite q-values,
#'   recompute q-values from \code{score_col} per run.
#' @param id_mode Character. \code{"id"} (default) uses \code{id = EG.PrecursorId} and keeps
#'   \code{run} separate; \code{"runxid"} uses \code{id = "run||precursor"} for global uniqueness.
#' @param unit Character. Unit metadata stored in the returned object.
#' @param scope Character. Scope metadata stored in the returned object.
#' @param q_source Character. Source label stored in metadata.
#'
#' @return
#' A \code{dfdr_tbl} (tibble subclass) with one row per \code{run} \eqn{\times}
#' \code{precursor}. The returned table contains the standard columns required by
#' \code{diagFDR}: \code{id}, \code{is_decoy}, \code{q}, \code{pep}, \code{run}, \code{score}.
#' Metadata (including whether q-values were recomputed) are stored in \code{attr(x, "meta")}.
#'
#' @details
#' When Spectronaut exports decoys with \code{NaN} q-values, recomputing q-values
#' from scores is performed *per run* to respect the run-wise hypothesis structure
#' and to avoid pooling scores across runs that may not be comparable.
#'
#' @examples
#' library(tibble)
#'
#' # Minimal toy Spectronaut-like report
#' rep <- tibble(
#'   R.FileName = c("run1","run1","run1","run2","run2","run2"),
#'   EG.PrecursorId = c("P1","P2","P3","P1","P2","P4"),
#'   EG.IsDecoy = c("False","True","False","False","True","False"),
#'   EG.Qvalue = c(0.01, NaN, 0.02, 0.03, NaN, 0.04),
#'   EG.PEP = c(0.02, 0.9, 0.05, 0.06, 0.95, 0.08),
#'   EG.Cscore = c(120, 10, 80, 100, 5, 60)
#' )
#'
#' x <- spectronaut_runxprecursor(
#'   rep,
#'   q_col = "EG.Qvalue",
#'   run_col = "R.FileName",
#'   recompute_q = TRUE,
#'   id_mode = "runxid"
#' )
#' x
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


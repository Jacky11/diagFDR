#' Internal: robust coercion of mokapot \code{is_target} column to logical
#'
#' @param x Vector containing mokapot \code{is_target} values (logical, numeric,
#'   or character).
#'
#' @return A logical vector of the same length as \code{x}.
#'
#' @keywords internal
mokapot_is_target_to_logical <- function(x) {
  # mokapot writes is_target as TRUE/FALSE (logical) or strings
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x == 1)
  z <- tolower(as.character(x))
  z %in% c("true", "t", "1", "yes", "y")
}

#' Compete winners by maximum score within keys (one per run+spectrum_id)
#'
#' @param df A data frame.
#' @param keys Character vector of column names defining groups to compete within.
#' @param score_col Character scalar giving the score column name.
#'
#' @return A data frame with at most one row per unique combination of \code{keys},
#'   retaining the row with the maximum score (ties broken by taking a single row).
#'
#' @keywords internal
compete_max_score <- function(df, keys, score_col) {
  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::slice_max(.data[[score_col]], n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
}

#' Read mokapot PSM text outputs (targets + decoys) and combine
#'
#' Reads the mokapot PSM output tables for targets and decoys (typically
#' \code{*.mokapot.psms.txt} and \code{*.mokapot.decoy.psms.txt}) and concatenates
#' them into a single table.
#'
#' @param target_path Character scalar. Path to the mokapot target PSM table
#'   (e.g. \code{*.mokapot.psms.txt}).
#' @param decoy_path Character scalar. Path to the mokapot decoy PSM table
#'   (e.g. \code{*.mokapot.decoy.psms.txt}).
#'
#' @return
#' A \link[tibble:tibble]{tibble} containing all rows from the target and decoy
#' tables. Column names and content are determined by mokapot.
#'
#' @examples
#' if (requireNamespace("readr", quietly = TRUE)) {
#'   # Create tiny example mokapot-like target/decoy tables and read them back
#'   tf1 <- tempfile(fileext = ".txt")
#'   tf2 <- tempfile(fileext = ".txt")
#'   txt1 <- "run\tspectrum_id\tmokapot score\tmokapot q-value\tmokapot PEP
#'   \tis_target\nr1\t1\t2.0\t0.01\t0.02\tTRUE\n"
#'   txt2 <- "run\tspectrum_id\tmokapot score\tmokapot q-value\tmokapot PEP
#'   \tis_target\nr1\t1\t1.0\t0.50\t0.60\tFALSE\n"
#'   writeLines(txt1, tf1)
#'   writeLines(txt2, tf2)
#'
#'   raw <- read_mokapot_psms(tf1, tf2)
#'   raw
#' }
#'
#' @export
read_mokapot_psms <- function(target_path, decoy_path) {
  if (!file.exists(target_path)) rlang::abort(paste0("File not found: ", target_path))
  if (!file.exists(decoy_path)) rlang::abort(paste0("File not found: ", decoy_path))

  T <- readr::read_tsv(target_path, show_col_types = FALSE)
  D <- readr::read_tsv(decoy_path, show_col_types = FALSE)

  dplyr::bind_rows(T, D)
}

#' Create a competed-winner universe from mokapot PSM outputs
#'
#' Constructs a "postprocessor universe" by selecting (competing) the single best
#' scoring entry (target or decoy) per \code{run + spectrum_id} from mokapot PSM
#' outputs. This yields one row per spectrum per run and can be used for downstream
#' target-decoy diagnostics.
#'
#' @param raw A combined mokapot table, typically the output of
#'   \code{\link{read_mokapot_psms}}.
#' @param run_col Character. Column name for run (default \code{"run"}).
#' @param spectrum_id_col Character. Column name for spectrum identifier (default
#'   \code{"spectrum_id"}).
#' @param score_col Character. Column name for mokapot score used for competition
#'   (default \code{"mokapot score"}; higher is better).
#' @param q_col Character. Column name for mokapot q-value (default
#'   \code{"mokapot q-value"}).
#' @param pep_col Character. Column name for mokapot PEP (default
#'   \code{"mokapot PEP"}).
#' @param id_mode Character. If \code{"runxid"} (default), sets \code{id} to
#'   \code{"run||spectrum_id"}. If \code{"id"}, sets \code{id = spectrum_id}.
#' @param unit Character. Unit metadata stored in the returned object.
#' @param scope Character. Scope metadata stored in the returned object.
#' @param q_source Character. Source label stored in metadata.
#' @param q_max_export Numeric. Export ceiling for q-values (typically 1 for mokapot).
#'
#' @return
#' An object of class \code{dfdr_tbl} (a tibble subclass) with one row per competed
#' \code{run + spectrum_id}. The returned table contains the standard columns
#' required by \code{diagFDR}: \code{id}, \code{is_decoy}, \code{q}, \code{pep},
#' \code{run}, and \code{score}. Metadata are stored in \code{attr(x, "meta")}.
#'
#' @examples
#' library(tibble)
#'
#' raw <- tibble(
#'   run = c("r1", "r1"),
#'   spectrum_id = c("1001", "1001"),
#'   `mokapot score` = c(2.0, 1.0),      # target wins
#'   `mokapot q-value` = c(0.01, 0.5),
#'   `mokapot PEP` = c(0.02, 0.6),
#'   is_target = c(TRUE, FALSE)
#' )
#'
#' x <- mokapot_competed_universe(raw, id_mode = "runxid")
#' x
#'
#' @export
mokapot_competed_universe <- function(raw,
                                      run_col = "run",
                                      spectrum_id_col = "spectrum_id",
                                      score_col = "mokapot score",
                                      q_col = "mokapot q-value",
                                      pep_col = "mokapot PEP",
                                      id_mode = c("runxid", "id"),
                                      unit = "psm",
                                      scope = "global",
                                      q_source = "mokapot q-value",
                                      q_max_export = 1) {
  id_mode <- match.arg(id_mode)

  req <- c(run_col, spectrum_id_col, score_col, q_col, pep_col, "is_target")
  miss <- setdiff(req, names(raw))
  if (length(miss)) rlang::abort(paste0("Missing required mokapot columns: ", paste(miss, collapse = ", ")))

  df <- raw |>
    dplyr::transmute(
      run = as.character(.data[[run_col]]),
      spectrum_id = as.character(.data[[spectrum_id_col]]),
      score = suppressWarnings(as.numeric(.data[[score_col]])),
      q = suppressWarnings(as.numeric(.data[[q_col]])),
      pep = suppressWarnings(as.numeric(.data[[pep_col]])),
      is_target = mokapot_is_target_to_logical(.data[["is_target"]])
    ) |>
    dplyr::mutate(
      is_decoy = !is_target
    ) |>
    dplyr::filter(
      !is.na(run), !is.na(spectrum_id),
      is.finite(score),
      is.finite(q), q >= 0, q <= 1,
      is.finite(pep), pep >= 0, pep <= 1
    )

  # compete winners per run+spectrum_id
  comp <- compete_max_score(df, keys = c("run", "spectrum_id"), score_col = "score")

  out <- comp |>
    dplyr::mutate(
      id = if (id_mode == "runxid") paste(run, spectrum_id, sep = "||") else spectrum_id
    ) |>
    dplyr::select(id, is_decoy, q, pep, run, score)

  as_dfdr_tbl(
    out,
    unit = unit,
    scope = scope,
    q_source = q_source,
    q_max_export = q_max_export,
    provenance = list(tool = "mokapot", id_mode = id_mode)
  )
}

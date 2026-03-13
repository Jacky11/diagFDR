#' Read mokapot PSM text outputs (targets + decoys) and combine
#'
#' @param target_path path to *.mokapot.psms.txt (targets)
#' @param decoy_path path to *.mokapot.decoy.psms.txt (decoys)
#' @return tibble combined
#' @export
read_mokapot_psms <- function(target_path, decoy_path) {
  if (!file.exists(target_path)) rlang::abort(paste0("File not found: ", target_path))
  if (!file.exists(decoy_path)) rlang::abort(paste0("File not found: ", decoy_path))

  T <- readr::read_tsv(target_path, show_col_types = FALSE)
  D <- readr::read_tsv(decoy_path, show_col_types = FALSE)

  dplyr::bind_rows(T, D)
}

#' Internal: robust coercion of mokapot is_target column to logical
#' @keywords internal
mokapot_is_target_to_logical <- function(x) {
  # mokapot writes is_target as TRUE/FALSE (logical) or strings
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x == 1)
  z <- tolower(as.character(x))
  z %in% c("true", "t", "1", "yes", "y")
}

#' Compete winners by max score within keys (one per run+spectrum_id)
#' @keywords internal
compete_max_score <- function(df, keys, score_col) {
  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::slice_max(.data[[score_col]], n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
}

#' Create a competed-winner universe from mokapot PSM outputs
#'
#' This is the generic "postprocessor" universe: one entry per run+spectrum_id
#' selected by the maximum mokapot score (target or decoy).
#'
#' @param raw combined mokapot table (from read_mokapot_psms)
#' @param run_col column name for run (default "run")
#' @param spectrum_id_col column name for spectrum_id (default "spectrum_id")
#' @param score_col column name for mokapot score (default "mokapot score")
#' @param q_col column name for mokapot q-value (default "mokapot q-value")
#' @param pep_col column name for mokapot PEP (default "mokapot PEP")
#' @param id_mode "runxid" makes id=run||spectrum_id (recommended); "id" makes id=spectrum_id
#' @param unit Unit metadata stored in the returned object.
#' @param scope Scope metadata stored in the returned object.
#' @param q_source Source label stored in metadata.
#' @param q_max_export Export ceiling for q-values (usually 1 for mokapot).
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

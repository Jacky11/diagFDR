## =====================================================================
##  Target/decoy label propagation when exported rows collapse to units
##
##  Implements the four rules of Supplementary Note S9 of the diagFDR
##  paper: any_decoy (decoy-if-any), winning_row, majority and split.
## =====================================================================

#' Internal: coerce a target/decoy column to logical
#'
#' Accepts logical, 0/1 numeric, integer64-like numeric, character
#' ("0"/"1"/"TRUE"/"FALSE") and factor. Factors are converted through
#' their labels, never through their level codes.
#'
#' @param v Vector.
#' @param where Column name used in error messages.
#' @return Logical vector (TRUE = decoy).
#' @keywords internal
as_decoy_logical <- function(v, where = "Decoy") {
  if (is.logical(v)) return(v)
  if (is.factor(v)) v <- as.character(v)
  if (is.character(v)) {
    u <- unique(v[!is.na(v)])
    ok <- c("0", "1", "TRUE", "FALSE", "True", "False", "true", "false")
    if (!all(u %in% ok))
      rlang::abort(paste0(where, ": unexpected character values: ",
                          paste(utils::head(u, 5), collapse = ", ")))
    return(ifelse(is.na(v), NA, v %in% c("1", "TRUE", "True", "true")))
  }
  vn <- suppressWarnings(as.numeric(v))
  u <- unique(vn[!is.na(vn)])
  if (!all(u %in% c(0, 1)))
    rlang::abort(paste0(where, ": expected 0/1, found ",
                        paste(utils::head(sort(u), 5), collapse = ", ")))
  vn == 1
}

#' Collapse exported rows to reported units under an explicit label-propagation rule
#'
#' Search engines export several rows per reported unit (one row per run and
#' precursor; many precursor rows per protein group). Building the list a
#' paper reports requires collapsing rows to units, and when the rows of a
#' unit do not all carry the same target/decoy label, a rule must decide
#' which label the unit inherits. This function implements the four rules
#' compared in the diagFDR paper (Supplementary Note S9):
#'
#' \describe{
#'   \item{\code{"split"}}{No label is propagated: the target rows and the
#'     decoy rows of an identifier are two distinct hypotheses, each kept at
#'     its own minimum q-value. A label-mixed identifier therefore yields one
#'     target unit \emph{and} one decoy unit. This reproduces DIA-NN's
#'     protein-group q-values and is the recommended rule wherever the
#'     software competes target and decoy versions of a group.}
#'   \item{\code{"winning_row"}}{The unit inherits the label of the row that
#'     supplied its minimum q-value. The defensible default when a genuine
#'     collapse to one unit per identifier is required.}
#'   \item{\code{"any_decoy"}}{The unit is a decoy if any of its rows is.
#'     What a direct \code{any(is_decoy)} aggregation produces; not
#'     recommended on label-mixed identifiers.}
#'   \item{\code{"majority"}}{The unit is a decoy if more than half of its
#'     rows are (ties resolve to target). Not recommended: its value depends
#'     on how many low-confidence rows survived the export filter.}
#' }
#'
#' In all rules the retained q-value (and PEP, when present) is the minimum
#' over the rows of the unit, and the retained score the maximum.
#'
#' @param rep A data.frame or tibble of exported rows (e.g. from
#'   \code{\link{read_diann_parquet}}).
#' @param key_col Name of the column identifying the reported unit
#'   (e.g. \code{"Precursor.Id"} or \code{"Protein.Group"}).
#' @param q_col Name of the q-value column to retain
#'   (e.g. \code{"Global.PG.Q.Value"}).
#' @param rule One of \code{"split"}, \code{"winning_row"},
#'   \code{"any_decoy"}, \code{"majority"}.
#' @param decoy_col Name of the target/decoy column (1/TRUE = decoy).
#' @param pep_col Name of the PEP column, or \code{NULL}. Ignored when
#'   absent from \code{rep}.
#' @param score_col Name of a score column, or \code{NULL}.
#' @param q_max_export Export ceiling of \code{q_col}, if known.
#' @param unit,scope,q_source Metadata stored in the result
#'   (see \code{\link{as_dfdr_tbl}}).
#' @param two_pass Logical: is this a two-pass (match-between-runs) export?
#' @param tool Name of the exporting tool, stored in provenance.
#' @param drop_empty_key Drop rows whose key is the empty string.
#'
#' @return A \code{dfdr_tbl} with one row per unit (two per label-mixed
#'   identifier under \code{"split"}), whose provenance records
#'   \code{label_rule}, \code{key}, the number of label-mixed identifiers
#'   (\code{n_mixed}) and the number of rows dropped.
#'
#' @seealso \code{\link{label_is_intrinsic}}, \code{\link{dfdr_label_rule_audit}},
#'   \code{\link{diann_global_proteingroup}}, \code{\link{dfdr_picked_protein}}.
#'
#' @examples
#' rep <- data.frame(
#'   Protein.Group = c("P1", "P1", "P1", "P1", "P2", "P3"),
#'   Decoy = c(1, 0, 0, 0, 0, 1),
#'   Global.PG.Q.Value = c(0.004, 0.02, 0.30, 0.60, 0.001, 0.05))
#' collapse_universe(rep, "Protein.Group", "Global.PG.Q.Value", rule = "split")
#' collapse_universe(rep, "Protein.Group", "Global.PG.Q.Value", rule = "winning_row")
#' @export
collapse_universe <- function(rep, key_col, q_col,
                              rule = c("split", "winning_row", "any_decoy", "majority"),
                              decoy_col = "Decoy", pep_col = "PEP", score_col = NULL,
                              q_max_export = NA_real_,
                              unit = "protein_group", scope = "global",
                              q_source = q_col, two_pass = NA,
                              tool = NA_character_, drop_empty_key = TRUE) {
  rule <- match.arg(rule)
  need <- c(key_col, decoy_col, q_col)
  miss <- setdiff(need, names(rep))
  if (length(miss))
    rlang::abort(paste0("collapse_universe(): missing columns ",
                        paste(miss, collapse = ", ")))

  key <- as.character(rep[[key_col]])
  q   <- suppressWarnings(as.numeric(rep[[q_col]]))
  dc  <- as_decoy_logical(rep[[decoy_col]], decoy_col)
  pep <- if (!is.null(pep_col) && pep_col %in% names(rep))
    suppressWarnings(as.numeric(rep[[pep_col]])) else NULL
  sc  <- if (!is.null(score_col) && score_col %in% names(rep))
    suppressWarnings(as.numeric(rep[[score_col]])) else NULL

  keep <- !is.na(key) & !is.na(dc) & is.finite(q) & q >= 0 & q <= 1
  if (drop_empty_key) keep <- keep & nzchar(key)
  n_drop <- sum(!keep)
  if (n_drop > 0 && n_drop / length(keep) > 0.5)
    rlang::warn(sprintf("collapse_universe(%s, %s): dropping %.0f%% of rows",
                        key_col, q_col, 100 * n_drop / length(keep)))
  key <- key[keep]; q <- q[keep]; dc <- dc[keep]
  if (!is.null(pep)) pep <- pep[keep]
  if (!is.null(sc))  sc  <- sc[keep]
  if (!length(key))
    rlang::abort(paste0("collapse_universe(", key_col, ", ", q_col,
                        "): no usable rows after filtering."))

  ## number of label-mixed identifiers (same for all rules)
  n_lab   <- tapply(dc, key, function(z) length(unique(z)))
  n_mixed <- sum(n_lab > 1L)

  grp_key <- if (rule == "split") paste0(ifelse(dc, "DECOY:", ""), key) else key
  o <- order(grp_key, q, method = "radix")
  grp_key <- grp_key[o]; q <- q[o]; dc <- dc[o]; key <- key[o]
  if (!is.null(pep)) pep <- pep[o]
  if (!is.null(sc))  sc  <- sc[o]

  first <- !duplicated(grp_key)
  idk   <- key[first]
  grp   <- cumsum(first)
  n_tot <- tabulate(grp, nbins = sum(first))

  lab <- if (rule %in% c("split", "winning_row")) {
    dc[first]
  } else {
    n_dec <- as.vector(rowsum(as.integer(dc), grp, reorder = FALSE))
    if (rule == "any_decoy") n_dec > 0L else (n_dec / n_tot) > 0.5
  }
  if (anyNA(lab))
    rlang::abort("collapse_universe(): NA label after filtering; please report.")

  pep_u <- if (!is.null(pep)) as.vector(tapply(pep, grp, safe_min)) else NA_real_
  sc_u  <- if (!is.null(sc))  as.vector(tapply(sc,  grp, safe_max)) else NA_real_

  out <- data.frame(id = idk, is_decoy = as.logical(lab), q = q[first],
                    pep = pep_u, run = NA_character_, score = sc_u,
                    stringsAsFactors = FALSE)

  as_dfdr_tbl(out, unit = unit, scope = scope, q_source = q_source,
              q_max_export = q_max_export,
              provenance = list(tool = tool, label_rule = rule, key = key_col,
                                n_mixed = n_mixed, n_rows_in = length(keep),
                                n_dropped = n_drop),
              two_pass = two_pass)
}

#' Is the target/decoy label intrinsic to the grouping key?
#'
#' If no identifier carries both a target and a decoy row, all four
#' label-propagation rules of \code{\link{collapse_universe}} coincide and the
#' rule does not need to be audited. On a DIA-NN export this is the case at
#' precursor level (\code{Precursor.Id}) but not at protein-group level
#' (\code{Protein.Group}).
#'
#' @inheritParams collapse_universe
#' @return A list with \code{intrinsic} (logical), \code{n_mixed} (number
#'   of label-mixed identifiers) and \code{n_keys}.
#' @examples
#' rep <- data.frame(Protein.Group = c("P1", "P1", "P2"), Decoy = c(1, 0, 0))
#' label_is_intrinsic(rep, "Protein.Group")
#' @export
label_is_intrinsic <- function(rep, key_col, decoy_col = "Decoy") {
  key <- as.character(rep[[key_col]])
  dc  <- as.integer(as_decoy_logical(rep[[decoy_col]], decoy_col))
  ok  <- !is.na(key) & !is.na(dc) & nzchar(key)
  n_lab <- tapply(dc[ok], key[ok], function(z) length(unique(z)))
  list(intrinsic = sum(n_lab > 1L) == 0L,
       n_mixed = sum(n_lab > 1L), n_keys = length(n_lab))
}

#' Audit the dependence of a reported list on the label-propagation rule
#'
#' Collapses the same export under each of the four rules of
#' \code{\link{collapse_universe}} and reports, at the threshold
#' \code{alpha}, the accepted target and decoy counts, the recomputed decoy
#' ratio, the decoy-support index and the FDP-CV benchmark. The spread across
#' rules is the unstated degree of freedom the diagFDR paper asks authors to
#' report (main-text Table 3, Panel B).
#'
#' @inheritParams collapse_universe
#' @param alpha Threshold at which the list is reported.
#' @param rules Character vector of rules to compare.
#' @return A data.frame with one row per rule and columns \code{rule},
#'   \code{n_units}, \code{T_alpha}, \code{D_alpha}, \code{decoy_ratio},
#'   \code{support_index}, \code{CV_benchmark}, \code{n_mixed}.
#' @examples
#' rep <- data.frame(
#'   Protein.Group = rep(paste0("P", 1:200), each = 3),
#'   Decoy = rbinom(600, 1, 0.1),
#'   Global.PG.Q.Value = runif(600, 0, 0.05))
#' dfdr_label_rule_audit(rep, "Protein.Group", "Global.PG.Q.Value", alpha = 0.01)
#' @export
dfdr_label_rule_audit <- function(rep, key_col, q_col, alpha = 0.01,
                                  rules = c("split", "winning_row", "any_decoy", "majority"),
                                  decoy_col = "Decoy", pep_col = "PEP",
                                  q_max_export = NA_real_, two_pass = NA) {
  rows <- lapply(rules, function(r) {
    u <- collapse_universe(rep, key_col, q_col, rule = r, decoy_col = decoy_col,
                           pep_col = pep_col, q_max_export = q_max_export,
                           two_pass = two_pass)
    acc <- u$q <= alpha
    Ta <- sum(acc & !u$is_decoy); Da <- sum(acc & u$is_decoy)
    data.frame(rule = r, n_units = nrow(u), T_alpha = Ta, D_alpha = Da,
               decoy_ratio = if (Ta > 0) Da / Ta else NA_real_,
               support_index = if (Da > 0) Da^(-0.5) else NA_real_,
               CV_benchmark = if (Da > 0) sqrt(2) * Da^(-0.5) else NA_real_,
               n_mixed = attr(u, "meta")$provenance$n_mixed,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  attr(out, "alpha") <- alpha
  out
}

#' Build a DIA-NN protein-group universe with an explicit label rule
#'
#' DIA-NN writes the same \code{Protein.Group} accession string for the
#' target and the decoy precursors of a group, so a protein-group list needs
#' a label-propagation rule (see \code{\link{collapse_universe}}). The
#' default \code{"split"} reproduces DIA-NN's protein-group q-values on the
#' exports analysed in the diagFDR paper.
#'
#' @param rep A DIA-NN report (see \code{\link{read_diann_parquet}}).
#' @param q_col \code{"Global.PG.Q.Value"} (experiment-wide list) or
#'   \code{"PG.Q.Value"} (run-wise values; taking the minimum across runs is
#'   the \eqn{\min_r} aggregation the paper flags as a scope mismatch).
#' @param rule Label-propagation rule.
#' @param protein_col Name of the protein-group column.
#' @param q_max_export Export ceiling of \code{q_col}, if known.
#' @param two_pass Logical: is this a \code{--reanalyse} export?
#' @return A \code{dfdr_tbl} at protein-group level.
#' @examples
#' rep <- data.frame(Protein.Group = c("P1", "P1", "P2"), Decoy = c(1, 0, 0),
#'                   Global.PG.Q.Value = c(0.2, 0.001, 0.003))
#' diann_global_proteingroup(rep)
#' @export
diann_global_proteingroup <- function(rep, q_col = "Global.PG.Q.Value",
                                      rule = c("split", "winning_row", "any_decoy", "majority"),
                                      protein_col = "Protein.Group",
                                      q_max_export = NA_real_, two_pass = NA) {
  rule <- match.arg(rule)
  diann_require_cols(rep, c(protein_col, "Decoy", q_col))
  scope <- if (q_col == "Global.PG.Q.Value") "global" else "min_over_runs"
  collapse_universe(rep, protein_col, q_col, rule = rule, decoy_col = "Decoy",
                    pep_col = NULL, q_max_export = q_max_export,
                    unit = "protein_group", scope = scope, q_source = q_col,
                    two_pass = diann_note_two_pass(two_pass), tool = "DIA-NN")
}

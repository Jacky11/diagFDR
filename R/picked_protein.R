## =====================================================================
##  Picked-protein competition (Savitski et al. 2015; The et al. 2021)
## =====================================================================

#' Picked-protein competition on a target/decoy list
#'
#' A procedural comparator for the label-propagation rules of
#' \code{\link{collapse_universe}}. Instead of inheriting a label, the
#' target and decoy versions of the same identifier compete directly: the
#' better of the pair (smaller q-value, or larger score) is kept and the
#' other discarded, and q-values are recomputed by target--decoy counting
#' over the retained entries. Because every retained entry comes from a
#' matched pair, the null is explicit and no label has to be propagated.
#'
#' The competition should be run on the \emph{full} (unfiltered) universe:
#' the decoy tail needed to recompute q-values is absent from an export
#' filtered at 1\%.
#'
#' @param x A \code{dfdr_tbl} (e.g. from \code{\link{collapse_universe}}
#'   with \code{rule = "split"}) or a data.frame with columns \code{id},
#'   \code{is_decoy} and \code{q} (and optionally \code{score}).
#' @param alpha Threshold at which the headline counts are reported.
#' @param by Evidence used for the competition and the recount:
#'   \code{"q"} (smaller wins) or \code{"score"} (larger wins).
#' @param add_decoy Added to the decoy count in the recomputed estimate:
#'   \code{0} gives \eqn{D/T} (as in the diagFDR paper), \code{1} the
#'   \eqn{(D+1)/T} convention.
#'
#' @return A list with \code{tbl} (a \code{dfdr_tbl} of the picked
#'   entries with the recomputed q-value in \code{q}, the input evidence in
#'   \code{q_in}), \code{T_alpha}, \code{D_alpha}, \code{decoy_ratio},
#'   \code{support_index}, \code{CV_benchmark}, \code{n_pairs} (identifiers
#'   carrying both labels), \code{n_ids}, and \code{accepted_targets}.
#'
#' @references Savitski MM et al. (2015) A scalable approach for protein false
#'   discovery rate estimation in large proteomic data sets. Mol Cell
#'   Proteomics 14:2394--2404.
#' @examples
#' df <- data.frame(id = c("P1", "P1", "P2", "P3", "P3"),
#'                  is_decoy = c(FALSE, TRUE, FALSE, TRUE, FALSE),
#'                  q = c(0.001, 0.2, 0.002, 0.004, 0.05))
#' dfdr_picked_protein(df, alpha = 0.01)
#' @export
dfdr_picked_protein <- function(x, alpha = 0.01, by = c("q", "score"), add_decoy = 0L) {
  by <- match.arg(by)
  df <- as.data.frame(x)
  if (!all(c("id", "is_decoy", "q") %in% names(df)))
    rlang::abort("`x` must have columns id, is_decoy, q.")
  if (by == "score" && !"score" %in% names(df))
    rlang::abort("by = \"score\" requires a `score` column.")
  ev <- if (by == "q") df$q else -df$score   # smaller is better in both cases
  ok <- !is.na(ev) & !is.na(df$is_decoy) & !is.na(df$id)
  df <- df[ok, , drop = FALSE]; ev <- ev[ok]

  ## best evidence per (id, label)
  grp  <- paste0(ifelse(df$is_decoy, "D:", "T:"), df$id)
  best <- tapply(ev, grp, min)
  bid  <- sub("^[DT]:", "", names(best))
  blab <- startsWith(names(best), "D:")
  qT <- stats::setNames(rep(Inf, length(unique(bid))), unique(bid))
  qD <- qT
  qT[bid[!blab]] <- best[!blab]
  qD[bid[blab]]  <- best[blab]
  ids <- names(qT)
  n_pairs <- sum(is.finite(qT) & is.finite(qD))

  pick_is_decoy <- qD < qT            # decoy survives only if strictly better
  pick_ev <- ifelse(pick_is_decoy, qD, qT)
  keep <- is.finite(pick_ev)
  P <- data.frame(id = ids[keep], is_decoy = unname(pick_is_decoy[keep]),
                  q_in = unname(pick_ev[keep]), stringsAsFactors = FALSE)
  if (by == "score") P$q_in <- -P$q_in

  ## recompute q-values by target--decoy counting; ties keep decoys first
  key <- if (by == "q") P$q_in else -P$q_in
  o  <- order(key, !P$is_decoy, method = "radix")
  P  <- P[o, , drop = FALSE]
  Tc <- cumsum(!P$is_decoy); Dc <- cumsum(P$is_decoy)
  fdr <- (Dc + add_decoy) / pmax(Tc, 1)
  P$q <- rev(cummin(rev(fdr)))
  P$pep <- NA_real_; P$run <- NA_character_
  P$score <- if (by == "score") P$q_in else NA_real_

  acc <- P$q <= alpha
  Ta <- sum(acc & !P$is_decoy); Da <- sum(acc & P$is_decoy)

  tbl <- as_dfdr_tbl(P[, c("id", "is_decoy", "q", "pep", "run", "score", "q_in")],
                     unit = "protein_group", scope = "picked",
                     q_source = paste0("picked_recount(", by, ")"),
                     q_max_export = 1,
                     provenance = list(method = "picked_protein", by = by,
                                       add_decoy = add_decoy, n_pairs = n_pairs))
  list(tbl = tbl, alpha = alpha, T_alpha = Ta, D_alpha = Da,
       decoy_ratio = if (Ta > 0) Da / Ta else NA_real_,
       support_index = if (Da > 0) Da^(-0.5) else NA_real_,
       CV_benchmark = if (Da > 0) sqrt(2) * Da^(-0.5) else NA_real_,
       n_pairs = n_pairs, n_ids = length(ids),
       n_pick_T = sum(!P$is_decoy), n_pick_D = sum(P$is_decoy),
       accepted_targets = P$id[acc & !P$is_decoy])
}

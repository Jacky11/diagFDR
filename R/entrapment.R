## =====================================================================
##  External validation by entrapment (A3)
##
##  Raw and size-corrected entrapment FDP estimates, Wilson intervals,
##  the scope-discordant-set error rate, entrapment attribution for
##  DIA-NN protein groups, and the entrapment-to-target size ratio r.
## =====================================================================

#' Entrapment-based FDP estimate of an accepted list
#'
#' Entrapment augments the target space with a proteome known to be absent
#' from the sample, searched as ordinary targets. Accepted matches mapping
#' only to it are empirical false positives independent of the decoys. With
#' \eqn{E_\alpha} such matches among the \eqn{T_\alpha} accepted targets
#' (\eqn{T_\alpha} counts every accepted non-decoy row, entrapment matches
#' included), the raw estimate \eqn{E_\alpha/T_\alpha} is the lower bound of
#' Wen et al. (2025) on the FDP of the combined target-plus-entrapment list;
#' their size-corrected ("combined") estimate \eqn{E_\alpha(1+1/r)/T_\alpha}
#' overestimates that FDP on average and is read as an upper bound, where
#' \eqn{r} is the entrapment-to-target size ratio measured at the level where
#' competition occurs (see \code{\link{dfdr_entrapment_size_ratio}}). The
#' "sample" estimate \eqn{(E_\alpha/r)/(T_\alpha-E_\alpha)} of the FDP among
#' the original targets only is also returned; Wen et al. show that it is
#' neither a lower nor an upper bound and typically underestimates.
#'
#' Wilson intervals are binomial-style descriptive intervals: they treat
#' accepted units as independent and do not account for clustering induced
#' by shared peptides, runs or protein groups.
#'
#' @param x A \code{dfdr_tbl}.
#' @param entrapment Either a logical vector of length \code{nrow(x)}
#'   (TRUE for identifications mapping only to the entrapment proteome) or
#'   a character vector of entrapment \code{id}s
#'   (see \code{\link{diann_entrapment_ids}}).
#' @param alpha Threshold(s) at which the list is reported.
#' @param r Entrapment-to-target size ratio, or \code{NA} to skip the
#'   corrected estimates.
#' @param conf Confidence level of the Wilson interval.
#'
#' @return A data.frame with one row per \code{alpha}: \code{T_alpha} (accepted
#'   targets, entrapment matches included), \code{E_alpha}, \code{fdp_raw},
#'   \code{lo}, \code{hi}, \code{r}, \code{fdp_corrected} (combined list),
#'   \code{fdp_corrected_target_only}, and \code{D_alpha}/\code{decoy_ratio}
#'   of the same list for comparison.
#'
#' @references Wen B, Freestone J, Riffle M, MacCoss MJ, Noble WS, Keich U (2025)
#'   Assessment of false discovery rate control in tandem mass spectrometry
#'   analysis using entrapment. Nat Methods.
#' @examples
#' sim <- dfdr_simulate_tdc(n = 5000, seed = 1)
#' ## pretend one target in ten is an entrapment sequence
#' ent <- !sim$is_decoy & (seq_len(nrow(sim)) %% 10 == 0) & attr(sim, "truth")$is_null
#' dfdr_entrapment(sim, ent, alpha = c(0.01, 0.05), r = 0.1)
#' @export
dfdr_entrapment <- function(x, entrapment, alpha = 0.01, r = NA_real_, conf = 0.95) {
  if (!inherits(x, "dfdr_tbl")) rlang::abort("`x` must be a dfdr_tbl.")
  is_ent <- .entrapment_flag(x, entrapment)
  out <- lapply(alpha, function(a) {
    acc <- !is.na(x$q) & x$q <= a
    Ta <- sum(acc & !x$is_decoy)
    Da <- sum(acc & x$is_decoy)
    Ea <- sum(acc & !x$is_decoy & is_ent)
    ci <- wilson_ci(Ea, Ta, conf = conf)
    data.frame(alpha = a, T_alpha = Ta, D_alpha = Da,
               decoy_ratio = if (Ta > 0) Da / Ta else NA_real_,
               E_alpha = Ea,
               fdp_raw = if (Ta > 0) Ea / Ta else NA_real_,
               lo = ci[1], hi = ci[2], r = r,
               fdp_corrected = if (is.finite(r) && Ta > 0)
                 Ea * (1 + 1 / r) / Ta else NA_real_,
               fdp_corrected_target_only = if (is.finite(r) && Ta - Ea > 0)
                 (Ea / r) / (Ta - Ea) else NA_real_)
  })
  do.call(rbind, out)
}

#' Internal: resolve the `entrapment` argument to a logical flag per row
#' @keywords internal
.entrapment_flag <- function(x, entrapment) {
  if (is.logical(entrapment)) {
    if (length(entrapment) != nrow(x))
      rlang::abort("logical `entrapment` must have length nrow(x).")
    entrapment[is.na(entrapment)] <- FALSE
    return(entrapment)
  }
  if (is.character(entrapment)) {
    ids <- sub("^.*\\|\\|", "", x$id)   # strip a "run||" prefix if present
    return(ids %in% entrapment | x$id %in% entrapment)
  }
  rlang::abort("`entrapment` must be a logical vector or a character vector of ids.")
}

#' Entrapment FDP on the discordant and shared parts of two accepted lists
#'
#' For two accepted target lists \eqn{A_1} and \eqn{A_2} on a common
#' identifier, computes the raw entrapment FDP on the shared core
#' \eqn{A_1 \cap A_2}, on the identifications added by \eqn{A_2}
#' (\eqn{A_2 \setminus A_1}), on those removed (\eqn{A_1 \setminus A_2}),
#' and on the symmetric difference (the \emph{scope-discordant-set error
#' rate} of the diagFDR paper, Equation 5). When neither list contains the
#' other, the added and removed sets should be read separately, since a
#' single rate can average an error-rich addition against an error-poor
#' removal.
#'
#' @param accepted1,accepted2 Character vectors of accepted target ids
#'   (see \code{\link{dfdr_accepted}}).
#' @param entrapment_ids Character vector of ids known to map only to the
#'   entrapment proteome.
#' @param conf Confidence level of the Wilson intervals.
#' @return A data.frame with rows \code{core}, \code{added}, \code{removed},
#'   \code{discordant}, \code{list1}, \code{list2}, columns \code{n},
#'   \code{E}, \code{fdp}, \code{lo}, \code{hi}; the attribute
#'   \code{"ratio"} holds \code{fdp[discordant] / fdp[core]} and
#'   \code{"jaccard"} the overlap of the two lists.
#' @examples
#' a1 <- paste0("p", 1:100); a2 <- paste0("p", 1:120)
#' dfdr_entrapment_discordant(a1, a2, entrapment_ids = c("p3", "p101", "p115"))
#' @export
dfdr_entrapment_discordant <- function(accepted1, accepted2, entrapment_ids, conf = 0.95) {
  a1 <- unique(as.character(accepted1)); a2 <- unique(as.character(accepted2))
  sets <- list(core = intersect(a1, a2), added = setdiff(a2, a1),
               removed = setdiff(a1, a2),
               discordant = union(setdiff(a1, a2), setdiff(a2, a1)),
               list1 = a1, list2 = a2)
  rows <- lapply(names(sets), function(nm) {
    s <- sets[[nm]]; n <- length(s); E <- sum(s %in% entrapment_ids)
    ci <- wilson_ci(E, n, conf = conf)
    data.frame(set = nm, n = n, E = E, fdp = if (n > 0) E / n else NA_real_,
               lo = ci[1], hi = ci[2], stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  attr(out, "ratio") <- out$fdp[out$set == "discordant"] / out$fdp[out$set == "core"]
  attr(out, "jaccard") <- jaccard_vec(a1, a2)
  attr(out, "nested") <- length(sets$removed) == 0 || length(sets$added) == 0
  out
}

#' Accepted target identifiers of a list at a threshold
#'
#' @param x A \code{dfdr_tbl}.
#' @param alpha Threshold.
#' @param strip_run Remove a \code{"run||"} prefix from ids (run x precursor
#'   universes) so that lists of different granularity index the same
#'   identifiers.
#' @return Character vector of accepted target ids.
#' @examples
#' sim <- dfdr_simulate_tdc(n = 2000, seed = 2)
#' length(dfdr_accepted(sim, 0.01))
#' @export
dfdr_accepted <- function(x, alpha = 0.01, strip_run = TRUE) {
  ids <- x$id[!x$is_decoy & !is.na(x$q) & x$q <= alpha]
  if (strip_run) ids <- sub("^.*\\|\\|", "", ids)
  unique(ids)
}

#' Entrapment status of protein-group strings
#'
#' Given accession strings such as DIA-NN's \code{Protein.Group}
#' (\code{"P12345;Q67890"}) and the set of entrapment accessions, flags each
#' string as an entrapment identification under the \code{"pure"} rule
#' (every accession belongs to the entrapment proteome; the conservative
#' rule used in the diagFDR paper) or the \code{"any"} rule (at least one
#' does).
#'
#' DIA-NN writes bare UniProt accessions with no species tag, so matching
#' must be done against the accession set and not a string prefix such as
#' \code{"ARATH"}; see \code{\link{read_fasta_accessions}}.
#'
#' @param groups Character vector of protein-group strings.
#' @param accessions Character vector of entrapment accessions.
#' @param rule \code{"pure"} or \code{"any"}.
#' @param sep Separator between accessions.
#' @return Logical vector.
#' @examples
#' dfdr_is_entrapment_group(c("A1;A2", "A1;P1", "P1"), accessions = c("A1", "A2"))
#' @export
dfdr_is_entrapment_group <- function(groups, accessions, rule = c("pure", "any"), sep = ";") {
  rule <- match.arg(rule)
  groups <- ifelse(is.na(groups), "", as.character(groups))
  parts <- strsplit(groups, sep, fixed = TRUE)
  vapply(parts, function(z) {
    z <- z[nzchar(z)]
    if (!length(z)) return(FALSE)
    if (rule == "pure") all(z %in% accessions) else any(z %in% accessions)
  }, logical(1))
}

#' Entrapment identifiers of a DIA-NN report
#'
#' Returns the unit identifiers (precursors by default, or protein groups)
#' whose protein group is an entrapment identification, for use as the
#' \code{entrapment} argument of \code{\link{dfdr_entrapment}} or the
#' \code{entrapment_ids} argument of \code{\link{dfdr_entrapment_discordant}}.
#'
#' @param rep A DIA-NN report.
#' @param accessions Character vector of entrapment accessions
#'   (\code{\link{read_fasta_accessions}}).
#' @param id_col Unit identifier column; use \code{"Protein.Group"} for
#'   protein-group universes.
#' @param protein_col Protein-group column.
#' @param rule \code{"pure"} or \code{"any"} (see
#'   \code{\link{dfdr_is_entrapment_group}}).
#' @return Character vector of ids. The function stops if no identifier is
#'   flagged, because an empty result is the silent-failure mode of a wrong
#'   accession set.
#' @examples
#' rep <- data.frame(Precursor.Id = c("AAK2", "BBR2", "CCK3"),
#'                   Protein.Group = c("A1", "A1;P1", "P2"))
#' diann_entrapment_ids(rep, accessions = "A1")
#' @export
diann_entrapment_ids <- function(rep, accessions, id_col = "Precursor.Id",
                                 protein_col = "Protein.Group", rule = c("pure", "any")) {
  rule <- match.arg(rule)
  diann_require_cols(rep, c(id_col, protein_col))
  flag <- dfdr_is_entrapment_group(rep[[protein_col]], accessions, rule = rule)
  ids <- unique(as.character(rep[[id_col]][flag]))
  if (!length(ids))
    rlang::abort(paste0("No identifier carries an entrapment accession. Check that ",
                        "`accessions` holds bare accessions (e.g. \"P12345\") ",
                        "matching the values of `", protein_col, "`."))
  ids
}

#' Accessions of a FASTA file
#'
#' Parses FASTA headers and returns the accession: the second field of
#' UniProt-style headers (\code{>sp|P12345|NAME_HUMAN ...}), otherwise the
#' first word after \code{>}.
#'
#' @param path Path to a FASTA file (plain or gzipped).
#' @param sequences Also return the sequences?
#' @return A character vector of accessions, or (with
#'   \code{sequences = TRUE}) a list with \code{ids} and \code{seqs}.
#' @examples
#' f <- tempfile(fileext = ".fasta")
#' writeLines(c(">sp|A0A1P8|X_ARATH test", "MKRPEPTIDEK", ">tr|B1|Y", "AAKR"), f)
#' read_fasta_accessions(f)
#' @export
read_fasta_accessions <- function(path, sequences = FALSE) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  hdr <- startsWith(lines, ">")
  heads <- lines[hdr]
  first <- sub("^>\\s*", "", heads)
  first <- sub("\\s.*$", "", first)
  parts <- strsplit(first, "|", fixed = TRUE)
  ids <- vapply(parts, function(z) if (length(z) >= 2) z[2] else z[1], character(1))
  if (!sequences) return(ids)
  grp <- cumsum(hdr)
  seqs <- vapply(split(lines[!hdr], grp[!hdr]), paste0, character(1), collapse = "")
  seqs <- toupper(gsub("\\s", "", seqs))
  list(ids = ids, seqs = unname(seqs[as.character(seq_along(ids))]))
}

#' Internal: in-silico digestion
#' @param seq Protein sequence.
#' @param sites Residues after which cleavage occurs.
#' @param no_cut_before Residues that block cleavage when following a site
#'   (\code{"P"} for the classical trypsin rule; \code{character(0)} for
#'   DIA-NN's \code{--cut K*,R*}).
#' @param missed Maximum number of missed cleavages.
#' @param min_len,max_len Peptide length range.
#' @keywords internal
.digest <- function(seq, sites = c("K", "R"), no_cut_before = character(0),
                    missed = 1L, min_len = 7L, max_len = 30L) {
  ch <- strsplit(seq, "")[[1]]
  n <- length(ch)
  if (!n) return(character(0))
  cut <- which(ch[-n] %in% sites)
  if (length(no_cut_before)) cut <- cut[!ch[cut + 1] %in% no_cut_before]
  bounds <- c(0L, cut, n)
  starts <- bounds[-length(bounds)] + 1L
  ends   <- bounds[-1]
  k <- length(starts)
  out <- character(0)
  for (m in 0:missed) {
    if (k - m < 1) break
    i <- seq_len(k - m)
    s <- starts[i]; e <- ends[i + m]
    len <- e - s + 1L
    keep <- len >= min_len & len <= max_len
    if (any(keep)) out <- c(out, substring(seq, s[keep], e[keep]))
  }
  unique(out)
}

#' Entrapment-to-target size ratio r
#'
#' Computes \eqn{r = N_{ET}/N_T} at protein level (number of entries) and at
#' peptide level (number of distinct in-silico tryptic peptides), with and
#' without the peptides shared between the two databases. The corrected
#' entrapment estimate of \code{\link{dfdr_entrapment}} depends on which
#' \eqn{r} is used, which is why both must be reported. Peptide-level
#' digestion of full proteomes takes a few minutes.
#'
#' @param target_fasta,entrapment_fasta Paths to the target (including
#'   contaminants, concatenated or as a character vector of paths) and
#'   entrapment FASTA files.
#' @param sites,no_cut_before,missed,min_len,max_len Digestion rule; the
#'   defaults reproduce DIA-NN \code{--cut K*,R* --missed-cleavages 1} with
#'   peptide length 7--30. Use \code{no_cut_before = "P"} for the classical
#'   trypsin rule.
#' @param peptide_level Also digest (slow on full proteomes)?
#' @return A data.frame with \code{quantity} and \code{value}: numbers of
#'   proteins and peptides, shared peptides, \code{r_protein},
#'   \code{r_peptide}, \code{r_peptide_unique}.
#' @examples
#' t <- tempfile(); e <- tempfile()
#' writeLines(c(">sp|P1|A", "MKAAAAAAAKRPEPTIDEKLLLLLLLR", ">sp|P2|B", "GGGGGGGGKRRRR"), t)
#' writeLines(c(">sp|A1|C", "MKAAAAAAAKRQQQQQQQQK"), e)
#' dfdr_entrapment_size_ratio(t, e)
#' @export
dfdr_entrapment_size_ratio <- function(target_fasta, entrapment_fasta,
                                       sites = c("K", "R"), no_cut_before = character(0),
                                       missed = 1L, min_len = 7L, max_len = 30L,
                                       peptide_level = TRUE) {
  tg <- lapply(target_fasta, read_fasta_accessions, sequences = TRUE)
  tg <- list(ids = unlist(lapply(tg, `[[`, "ids")), seqs = unlist(lapply(tg, `[[`, "seqs")))
  en <- read_fasta_accessions(entrapment_fasta, sequences = TRUE)
  n_t <- length(tg$ids); n_e <- length(en$ids)
  out <- data.frame(quantity = c("target_proteins", "entrapment_proteins", "r_protein"),
                    value = c(n_t, n_e, n_e / n_t), stringsAsFactors = FALSE)
  if (!peptide_level) return(out)
  dig <- function(s) unique(unlist(lapply(s, .digest, sites = sites,
                                          no_cut_before = no_cut_before, missed = missed,
                                          min_len = min_len, max_len = max_len)))
  pt <- dig(tg$seqs); pe <- dig(en$seqs)
  sh <- intersect(pt, pe)
  rbind(out, data.frame(
    quantity = c("target_peptides", "entrapment_peptides", "shared_peptides",
                 "r_peptide", "r_peptide_unique"),
    value = c(length(pt), length(pe), length(sh),
              length(pe) / length(pt), length(setdiff(pe, sh)) / length(pt)),
    stringsAsFactors = FALSE))
}

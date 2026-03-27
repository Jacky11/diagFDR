#' Read mzIdentML into a \code{dfdr_tbl} (generic; score-based TDC q-values)
#'
#' Extracts a competed PSM universe (rank-1 by default) from an mzIdentML file,
#' determines target/decoy labels, selects a single numeric PSM score CV term,
#' and reconstructs q-values using target-decoy counting (TDC).
#'
#' This function is intended for workflows where mzIdentML does not provide
#' explicit q-values or PEPs.
#'
#' @param mzid_path Character scalar. Path to an mzIdentML file (\code{.mzid}).
#' @param rank Integer scalar. Which \code{SpectrumIdentificationItem@rank} to use
#'   (default 1).
#' @param score_accession_preference Character vector. Preferred PSI-MS CV accessions
#'   to treat as the primary PSM score, in priority order.
#' @param score_direction One of \code{"auto"}, \code{"lower_better"}, \code{"higher_better"}.
#'   If \code{"auto"}, applies simple rules for common e-value/expectation accessions.
#' @param add_decoy Integer scalar. Additive correction in the TDC FDR estimate:
#'   \eqn{\widehat{FDR} = (D + add\_decoy) / T}. Default 1.
#' @param min_score_coverage Numeric in \eqn{(0,1]}. Minimum fraction of PSMs required to
#'   have the chosen score CV term. Default 1.0 (strict; no missing scores allowed).
#' @param decoy_regex Character scalar. Regex used to infer decoys from protein accessions
#'   if \code{PeptideEvidence@isDecoy} is not informative.
#' @param unit Character. Stored in \code{dfdr_tbl} metadata (default \code{"psm"}).
#' @param scope Character. Stored in \code{dfdr_tbl} metadata (default \code{NA}).
#' @param provenance Named list. Stored in metadata (tool, version, parameters, command, etc.).
#'
#' @return
#' A \code{dfdr_tbl} (tibble subclass) with one row per extracted PSM (at the requested
#' \code{rank}). The returned object contains:
#' \itemize{
#'   \item \code{id}: a PSM identifier derived from run and spectrum ID (\code{"run||spectrumID"});
#'   \item \code{is_decoy}: logical target/decoy label;
#'   \item \code{score}: an internal score where larger values mean better matches;
#'   \item \code{q}: reconstructed monotone q-values from TDC using \code{score} and \code{is_decoy};
#'   \item \code{pep}: \code{NA} (mzIdentML PEPs are not parsed here);
#'   \item \code{run}: run identifier when available.
#' }
#' Metadata are stored in \code{attr(x, "meta")} (including \code{unit}, \code{scope},
#' \code{q_source}, and \code{provenance}).
#'
#' @examples
#' if (requireNamespace("xml2", quietly = TRUE)) {
#'   # Minimal mzIdentML-like file sufficient for diagFDR's parser:
#'   # - 1 target and 1 decoy PSM at rank=1, with a numeric cvParam score
#'   tmp <- tempfile(fileext = ".mzid")
#'   mzid_txt <- paste0(
#'     "<?xml version='1.0' encoding='UTF-8'?>\n",
#'     "<MzIdentML xmlns='http://psidev.info/psi/pi/mzIdentML/1.1'>\n",
#'     "  <SequenceCollection>\n",
#'     "    <DBSequence id='DBSeq_t' accession='PROT1'/>\n",
#'     "    <DBSequence id='DBSeq_d' accession='REV_PROT2'/>\n",
#'     "    <PeptideEvidence id='PE_t' dBSequence_ref='DBSeq_t' isDecoy='false'/>\n",
#'     "    <PeptideEvidence id='PE_d' dBSequence_ref='DBSeq_d' isDecoy='true'/>\n",
#'     "  </SequenceCollection>\n",
#'     "  <DataCollection>\n",
#'     "    <AnalysisData>\n",
#'     "      <SpectrumIdentificationList>\n",
#'     "        <SpectrumIdentificationResult spectraData_ref='runA' spectrumID='scan=1'>\n",
#'     "          <SpectrumIdentificationItem rank='1'>\n",
#'     "            <PeptideEvidenceRef peptideEvidence_ref='PE_t'/>\n",
#'     "            <cvParam accession='MS:1001331' name='X!Tandem:hyperscore' value='50.0'/>\n",
#'     "          </SpectrumIdentificationItem>\n",
#'     "        </SpectrumIdentificationResult>\n",
#'     "        <SpectrumIdentificationResult spectraData_ref='runA' spectrumID='scan=2'>\n",
#'     "          <SpectrumIdentificationItem rank='1'>\n",
#'     "            <PeptideEvidenceRef peptideEvidence_ref='PE_d'/>\n",
#'     "            <cvParam accession='MS:1001331' name='X!Tandem:hyperscore' value='10.0'/>\n",
#'     "          </SpectrumIdentificationItem>\n",
#'     "        </SpectrumIdentificationResult>\n",
#'     "      </SpectrumIdentificationList>\n",
#'     "    </AnalysisData>\n",
#'     "  </DataCollection>\n",
#'     "</MzIdentML>\n"
#'   )
#'   writeLines(mzid_txt, tmp)
#'
#'   x <- read_dfdr_mzid(tmp, rank = 1L, score_direction = "higher_better")
#'   x
#'   range(x$q)
#' }
#'
#' @export
#'
#' @importFrom tibble tibble
#' @importFrom rlang abort
read_dfdr_mzid <- function(
    mzid_path,
    rank = 1L,
    score_accession_preference = c(
      # expectation/e-values first (lower is better)
      "MS:1002257", # Comet:expectation value
      "MS:1001330", # X!Tandem:expect
      "MS:1001328", # OMSSA:evalue
      "MS:1002052", # MS-GF:SpecEValue

      # then higher-better scores
      "MS:1002049", # MS-GF:RawScore
      "MS:1001331", # X!Tandem:hyperscore
      "MS:1001171", # Mascot:score
      "MS:1001950", # PEAKS:peptideScore

      # keep PeptideShaker score last (can saturate)
      "MS:1002466"  # PeptideShaker PSM score
    ),
    score_direction = c("auto", "lower_better", "higher_better"),
    add_decoy = 1L,
    min_score_coverage = 1.0,
    decoy_regex = "(^##|_REVERSED$|^REV_|^DECOY_)",
    unit = "psm",
    scope = NA_character_,
    provenance = list()
) {
  score_direction <- match.arg(score_direction)

  psms <- mzid_parse_psms(
    mzid_path = mzid_path,
    rank = rank,
    score_accession_preference = score_accession_preference,
    min_score_coverage = min_score_coverage,
    decoy_regex = decoy_regex
  )

  if (is.na(psms$score_accession[1])) {
    rlang::abort("No numeric score CV term found in mzIdentML; cannot compute q-values.")
  }

  score_dir <- mzid_score_direction(
    score_accession = psms$score_accession[1],
    score_direction = score_direction
  )

  score_internal <- mzid_to_internal_score(
    score_value = psms$score_value,
    score_direction = score_dir,
    score_accession = psms$score_accession[1]
  )

  # strict requirements for dfdr_tbl
  if (anyNA(psms$is_decoy)) rlang::abort("`is_decoy` contains NA after parsing; please fix upstream.")
  if (!any(psms$is_decoy)) rlang::abort("No decoys detected; cannot compute q-values.")
  if (!any(!psms$is_decoy)) rlang::abort("No targets detected; cannot compute q-values.")
  if (anyNA(score_internal)) rlang::abort("Chosen score is missing (NA) for some PSMs; cannot compute q-values.")
  if (all(is.na(score_internal))) rlang::abort("Score values are all NA; cannot compute q-values.")

  q <- tdc_qvalues(score = score_internal, is_decoy = psms$is_decoy, add_decoy = add_decoy)

  if (any(!is.finite(q))) rlang::abort("Reconstructed q contains non-finite values.")
  if (any(q < 0 | q > 1)) rlang::abort("Reconstructed q not in [0,1].")

  out <- tibble::tibble(
    id = psms$id,
    is_decoy = psms$is_decoy,
    q = q,
    pep = rep(NA_real_, length(q)),
    run = psms$run %||% rep(NA_character_, length(q)),
    score = score_internal
  )

  q_source <- paste0(
    "reconstructed_TDC_from_", psms$score_accession[1],
    " (", psms$score_name[1] %||% "unknown_score_name",
    "; add_decoy=", add_decoy, ")"
  )

  as_dfdr_tbl(
    out,
    unit = unit,
    scope = scope,
    q_source = q_source,
    q_max_export = NA_real_,
    provenance = provenance
  )
}


#' Parse rank-PSMs from mzIdentML (internal)
#'
#' Extracts SpectrumIdentificationItems at the requested rank, derives decoy labels
#' using PeptideEvidence@isDecoy when available, with fallback to DBSequence accession
#' regex inference, and selects a single numeric score term with sufficient coverage.
#'
#' @param mzid_path Character scalar. Path to an mzIdentML file.
#' @param rank Integer scalar. Which `SpectrumIdentificationItem@rank` to extract.
#' @param score_accession_preference Character vector of preferred score accessions.
#' @param min_score_coverage Numeric in (0,1]. Minimum fraction of PSMs required to have
#'   the chosen score CV term.
#' @param decoy_regex Character scalar. Regex for inferring decoys from DBSequence accessions.
#'
#' @return Tibble with columns id, run, is_decoy, score_accession, score_name, score_value.
#' @keywords internal
#'
#' @importFrom xml2 read_xml xml_find_all xml_parent xml_attr
#' @importFrom tibble tibble
#' @importFrom dplyr filter
#' @importFrom rlang abort
#' @return A \link[tibble:tibble]{tibble} with columns \code{id}, \code{run},
#'   \code{is_decoy}, \code{score_accession}, \code{score_name}, \code{score_value}.
#' @keywords internal
mzid_parse_psms <- function(
    mzid_path,
    rank = 1L,
    score_accession_preference = c("MS:1002257"),
    min_score_coverage = 1.0,
    decoy_regex = "(^##|_REVERSED$|^REV_|^DECOY_)"
) {
  if (!is.character(mzid_path) || length(mzid_path) != 1) {
    rlang::abort("`mzid_path` must be a character scalar.")
  }
  if (!file.exists(mzid_path)) rlang::abort("mzIdentML file not found: ", mzid_path)
  if (!is.numeric(min_score_coverage) || length(min_score_coverage) != 1 ||
      !is.finite(min_score_coverage) || min_score_coverage <= 0 || min_score_coverage > 1) {
    rlang::abort("`min_score_coverage` must be a single number in (0,1].")
  }

  doc <- xml2::read_xml(mzid_path)

  # ---- DBSequence lookup: id -> accession (for decoy inference fallback) ----
  dbseq_nodes <- xml2::xml_find_all(doc, ".//*[local-name()='DBSequence']")
  db_lut <- tibble::tibble(
    dbsequence_id = xml2::xml_attr(dbseq_nodes, "id"),
    accession = xml2::xml_attr(dbseq_nodes, "accession")
  )

  # ---- PeptideEvidence lookup: id -> isDecoy + dBSequence_ref ----
  pe_nodes <- xml2::xml_find_all(doc, ".//*[local-name()='PeptideEvidence']")
  if (length(pe_nodes) == 0) {
    rlang::abort("No <PeptideEvidence> nodes found; cannot determine targets/decoys.")
  }

  pe_lut <- tibble::tibble(
    peptideEvidence_ref = xml2::xml_attr(pe_nodes, "id"),
    is_decoy_attr = tolower(xml2::xml_attr(pe_nodes, "isDecoy")) %in% "true",
    dbsequence_ref = xml2::xml_attr(pe_nodes, "dBSequence_ref")
  )

  # join to get accession for each PeptideEvidence
  pe_lut$accession <- db_lut$accession[match(pe_lut$dbsequence_ref, db_lut$dbsequence_id)]

  # ---- Find rank-n SpectrumIdentificationItem nodes ----
  sii_xpath <- sprintf(".//*[local-name()='SpectrumIdentificationItem'][@rank='%s']",
                       as.integer(rank))
  sii_nodes <- xml2::xml_find_all(doc, sii_xpath)
  if (length(sii_nodes) == 0) {
    rlang::abort(paste0("No SpectrumIdentificationItem nodes found for rank=", rank))
  }

  sir_nodes <- xml2::xml_parent(sii_nodes)
  run <- xml2::xml_attr(sir_nodes, "spectraData_ref")
  spectrumID <- xml2::xml_attr(sir_nodes, "spectrumID")
  id <- paste0(run, "||", spectrumID)

  # ---- Evidence refs per SII ----
  pe_refs <- lapply(sii_nodes, function(n) {
    refs <- xml2::xml_find_all(n, ".//*[local-name()='PeptideEvidenceRef']")
    xml2::xml_attr(refs, "peptideEvidence_ref")
  })

  # ---- Decoy determination: attribute first, fallback to accession regex ----
  is_decoy_attr_based <- vapply(pe_refs, function(refs) {
    hit <- pe_lut$is_decoy_attr[match(refs, pe_lut$peptideEvidence_ref)]
    if (length(hit) == 0 || all(is.na(hit))) return(NA)
    any(hit %in% TRUE, na.rm = TRUE)
  }, logical(1))

  has_any_decoy_attr <- any(is_decoy_attr_based %in% TRUE, na.rm = TRUE)

  if (has_any_decoy_attr) {
    # Use isDecoy attribute (standard) if it yields at least one decoy
    is_decoy <- is_decoy_attr_based
  } else {
    # Fallback: infer decoys from DBSequence accession via regex
    is_decoy <- vapply(pe_refs, function(refs) {
      acc <- pe_lut$accession[match(refs, pe_lut$peptideEvidence_ref)]
      acc <- acc[!is.na(acc)]
      if (length(acc) == 0) {
        rlang::abort("Cannot infer decoy status: missing DBSequence accession for some PeptideEvidence.")
      }
      any(grepl(decoy_regex, acc, perl = TRUE))
    }, logical(1))
  }

  if (anyNA(is_decoy)) {
    rlang::abort("`is_decoy` contains NA after parsing/inference.")
  }

  # ---- Extract numeric cvParams per SII ----
  cv_df_list <- lapply(sii_nodes, function(n) {
    cvp <- xml2::xml_find_all(n, ".//*[local-name()='cvParam']")
    if (length(cvp) == 0) {
      return(tibble::tibble(accession = character(), name = character(), value = numeric()))
    }
    df <- tibble::tibble(
      accession = xml2::xml_attr(cvp, "accession"),
      name = xml2::xml_attr(cvp, "name"),
      value = suppressWarnings(as.numeric(xml2::xml_attr(cvp, "value")))
    )
    dplyr::filter(df, !is.na(.data$value))
  })

  # Compute coverage per accession (fraction of SIIs having it)
  acc_per_sii <- lapply(cv_df_list, function(df) unique(df$accession))
  all_acc <- sort(unique(unlist(acc_per_sii)))
  all_acc <- all_acc[!is.na(all_acc)]

  if (length(all_acc) == 0) {
    rlang::abort("No numeric cvParam values found on SpectrumIdentificationItem.")
  }

  coverage <- vapply(all_acc, function(a) {
    mean(vapply(acc_per_sii, function(v) a %in% v, logical(1)))
  }, numeric(1))

  eligible <- all_acc[coverage[all_acc] >= min_score_coverage]

  if (length(eligible) == 0) {
    # Strict behaviour: refuse to compute q if no single score term is widely present
    rlang::abort(
      "No single numeric score cvParam meets min_score_coverage=",
      min_score_coverage, ". ",
      "Cannot compute q-values without dropping rows or mixing score types."
    )
  }

  pref <- intersect(score_accession_preference, eligible)
  if (length(pref) > 0) {
    chosen_acc <- pref[1]
  } else {
    # Choose highest coverage among eligible; break ties by lexical order for determinism
    cov_elig <- coverage[eligible]
    chosen_acc <- eligible[order(cov_elig, decreasing = TRUE, eligible)][1]
  }

  # Pull that score for each SII; must not be missing
  score_value <- rep(NA_real_, length(sii_nodes))
  score_name  <- rep(NA_character_, length(sii_nodes))

  for (i in seq_along(cv_df_list)) {
    df <- cv_df_list[[i]]
    j <- which(df$accession == chosen_acc)
    if (length(j)) {
      score_value[i] <- df$value[j[1]]
      score_name[i]  <- df$name[j[1]]
    }
  }

  if (anyNA(score_value)) {
    rlang::abort(
      "Chosen score accession ", chosen_acc, " is missing for some PSMs (NA score_value). ",
      "Refusing to drop rows or mix score types. Consider lowering min_score_coverage."
    )
  }

  tibble::tibble(
    id = id,
    run = run,
    is_decoy = is_decoy,
    score_accession = chosen_acc,
    score_name = score_name,
    score_value = score_value
  )
}


#' Determine score direction (internal)
#'
#' @param score_accession Character. PSI-MS accession for the score term.
#' @param score_direction One of `"auto"`, `"lower_better"`, `"higher_better"`.
#'
#' @return `"lower_better"` or `"higher_better"`.
#' @keywords internal
mzid_score_direction <- function(score_accession,
                                 score_direction = c("auto", "lower_better", "higher_better")) {
  score_direction <- match.arg(score_direction)
  if (score_direction != "auto") return(score_direction)

  # E-values/expectation values are lower-better
  if (score_accession %in% c(
    "MS:1002257", # Comet:expectation value
    "MS:1001330", # X!Tandem:expect
    "MS:1001328", # OMSSA:evalue
    "MS:1002052"  # MS-GF:SpecEValue
  )) {
    return("lower_better")
  }

  "higher_better"
}


#' Convert raw score to an internal "higher is better" score (internal)
#'
#' @param score_value Numeric vector. Raw score values.
#' @param score_direction Character. `"lower_better"` or `"higher_better"`.
#' @param score_accession Character. Score accession (used for special handling).
#'
#' @return Numeric vector. Larger values mean better matches.
#' @keywords internal
mzid_to_internal_score <- function(score_value, score_direction, score_accession) {
  if (score_direction == "higher_better") return(score_value)

  # Use -log10 for common p/e-value like scores (strictly positive)
  if (score_accession %in% c("MS:1002257", "MS:1001330", "MS:1001328", "MS:1002052")) {
    eps <- .Machine$double.xmin
    return(-log10(pmax(score_value, eps)))
  }

  # Generic fallback
  -score_value
}


#' Compute monotone q-values by target-decoy counting (internal)
#'
#' Computes FDR estimates along the score-sorted list and returns monotone q-values:
#' \eqn{q(i) = \min_{j \ge i} \widehat{FDR}(j)}.
#'
#' @param score Numeric vector. Higher is better.
#' @param is_decoy Logical vector. TRUE=decoy, FALSE=target.
#' @param add_decoy Integer scalar. FDR hat = (D + add_decoy) / T.
#'
#' @return Numeric vector of q-values in `[0,1]`.
#' @keywords internal
tdc_qvalues <- function(score, is_decoy, add_decoy = 1L) {
  if (!is.numeric(score)) rlang::abort("`score` must be numeric.")
  if (!is.logical(is_decoy)) rlang::abort("`is_decoy` must be logical.")
  if (length(score) != length(is_decoy)) rlang::abort("`score` and `is_decoy` must have equal length.")
  if (anyNA(score)) rlang::abort("`score` contains NA; cannot compute q-values.")
  if (anyNA(is_decoy)) rlang::abort("`is_decoy` contains NA; cannot compute q-values.")

  ord <- order(score, decreasing = TRUE)

  dec <- is_decoy[ord]
  tar <- !is_decoy[ord]

  D <- cumsum(dec)
  T <- cumsum(tar)

  fdr_hat <- rep(NA_real_, length(ord))
  ok <- T > 0
  fdr_hat[ok] <- (D[ok] + add_decoy) / T[ok]

  q_ord <- rev(cummin(rev(fdr_hat)))

  q <- rep(NA_real_, length(score))
  q[ord] <- q_ord

  pmin(pmax(q, 0), 1)
}

#' Read MaxQuant \code{msms.txt} into a \code{dfdr_tbl} (PSM-level; reconstructed TDC q-values)
#'
#' Reads a MaxQuant \code{msms.txt} file and returns a \code{dfdr_tbl} at the PSM
#' level. q-values are reconstructed by target-decoy counting (TDC) using MaxQuant
#' \code{Score} (higher is better) and the \code{Reverse} column as the decoy
#' indicator (\code{"+"} for decoy).
#'
#' @details
#' Required columns in \code{msms.txt}:
#' \itemize{
#'   \item \code{id} (unique PSM identifier; numeric or character)
#'   \item \code{Raw file} (run name)
#'   \item \code{Reverse} (MaxQuant decoy indicator; \code{"+"} for decoy; blank/NA for target)
#'   \item \code{Score} (MaxQuant score; higher is better)
#'   \item \code{PEP} (posterior error probability; optional in practice, see \code{pep_mode})
#' }
#'
#' The function computes q-values by target-decoy counting (TDC) using \code{Score}
#' and the \code{Reverse} decoy label:
#' \deqn{\widehat{\mathrm{FDR}}(i) = (D(i) + add\_decoy) / T(i), \qquad
#'      q(i) = \min_{j \ge i}\widehat{\mathrm{FDR}}(j)}
#'
#' Contaminants are not treated as decoys. If the column
#' \code{Potential contaminant} exists and \code{exclude_contaminants = TRUE},
#' those rows are removed prior to computing q-values.
#'
#' @param path Character scalar. Path to MaxQuant \code{msms.txt}.
#' @param pep_mode Character. How to handle the \code{PEP} column:
#'   \describe{
#'     \item{\code{"drop"}}{Set \code{pep = NA} for all rows.}
#'     \item{\code{"sanitize"}}{Keep finite PEP in \code{[0,1]}; set other finite values to \code{NA} and warn.}
#'     \item{\code{"strict"}}{Error if any finite PEP is outside \code{[0,1]}.}
#'   }
#' @param exclude_contaminants Logical. If \code{TRUE} and a
#'   \code{Potential contaminant} column exists, rows with \code{"+"} are removed.
#'   Default is \code{TRUE}.
#' @param add_decoy Integer scalar. Additive correction in the TDC estimate:
#'   \eqn{\widehat{\mathrm{FDR}} = (D + add\_decoy)/T}. Default is 1.
#' @param unit Character. Unit stored in the returned object metadata (default \code{"psm"}).
#' @param scope Character. Scope stored in the returned object metadata (default \code{"global"}).
#' @param provenance Named list. Stored in the returned object metadata.
#'
#' @return
#' A \code{dfdr_tbl} (tibble subclass) with one row per PSM and columns:
#' \code{id}, \code{run}, \code{is_decoy}, \code{q}, \code{pep}, \code{score}.
#' The q-values are reconstructed by TDC from \code{Score} and \code{Reverse}.
#' Metadata are stored in \code{attr(x, "meta")} (including \code{unit}, \code{scope},
#' \code{q_source}, and \code{provenance}).
#'
#' @examples
#' if (requireNamespace("readr", quietly = TRUE)) {
#'   # Create a tiny MaxQuant-like msms.txt and read it
#'   tmp <- tempfile(fileext = ".txt")
#'   txt <- paste0(
#'     "id\tRaw file\tReverse\tScore\tPEP\tPotential contaminant\n",
#'     "1\trun1\t\t100\t0.01\t\n",
#'     "2\trun1\t+\t90\t0.80\t\n",
#'     "3\trun1\t\t80\t0.02\t+\n",   # contaminant row removed when exclude_contaminants=TRUE
#'     "4\trun2\t\t70\t0.05\t\n",
#'     "5\trun2\t+\t60\t0.90\t\n"
#'   )
#'   writeLines(txt, tmp)
#'
#'   x <- read_dfdr_maxquant_msms(tmp, pep_mode = "sanitize",
#'         exclude_contaminants = TRUE, add_decoy = 1L)
#'   x
#'   head(x$q)
#' }
#'
#' @export
read_dfdr_maxquant_msms <- function(path,
                                    pep_mode = c("drop", "sanitize", "strict"),
                                    exclude_contaminants = TRUE,
                                    add_decoy = 1L,
                                    unit = "psm",
                                    scope = "global",
                                    provenance = list()) {
  pep_mode <- match.arg(pep_mode)

  if (!is.character(path) || length(path) != 1) {
    rlang::abort("`path` must be a character scalar.")
  }
  if (!file.exists(path)) {
    rlang::abort(paste0("File not found: ", path))
  }

  x <- readr::read_tsv(
    file = path,
    show_col_types = FALSE,
    progress = FALSE,
    na = c("", "NA")
  )

  req <- c("id", "Raw file", "Reverse", "Score", "PEP")
  miss <- setdiff(req, names(x))
  if (length(miss)) {
    rlang::abort(paste0("msms.txt missing required columns: ", paste(miss, collapse = ", ")))
  }

  # Optional contaminant filtering
  if (isTRUE(exclude_contaminants) && ("Potential contaminant" %in% names(x))) {
    x <- dplyr::filter(x, is.na(.data[["Potential contaminant"]]) | .data[["Potential contaminant"]] != "+")
  }

  # Reverse -> is_decoy: "+" => TRUE; NA/blank => FALSE; anything else => error
  reverse_chr <- as.character(x[["Reverse"]])
  bad <- !is.na(reverse_chr) & nzchar(reverse_chr) & reverse_chr != "+"
  if (any(bad)) {
    vals <- sort(unique(reverse_chr[bad]))
    rlang::abort(paste0(
      "Unexpected values in `Reverse` column (expected '+' for decoy, blank/NA for target): ",
      paste(vals, collapse = ", ")
    ))
  }
  is_decoy <- ifelse(!is.na(reverse_chr) & reverse_chr == "+", TRUE, FALSE)

  # Score must be numeric + finite
  score_chr <- as.character(x[["Score"]])
  score <- suppressWarnings(as.numeric(score_chr))
  if (any(!is.finite(score))) {
    bad_sc <- which(!is.finite(score))
    ex <- unique(score_chr[bad_sc])
    ex <- ex[seq_len(min(8, length(ex)))]
    rlang::abort(paste0(
      "`Score` contains non-finite or non-numeric values; cannot compute q-values. Examples: ",
      paste(ex, collapse = ", ")
    ))
  }

  # PEP handling
  pep_chr <- as.character(x[["PEP"]])
  pep <- suppressWarnings(as.numeric(pep_chr))
  pep[!is.finite(pep)] <- NA_real_

  if (pep_mode == "drop") {
    pep <- rep(NA_real_, length(score))
  } else {
    bad_pep <- which(is.finite(pep) & (pep < 0 | pep > 1))
    if (length(bad_pep)) {
      if (pep_mode == "strict") {
        ex <- unique(pep_chr[bad_pep])
        ex <- ex[seq_len(min(8, length(ex)))]
        rlang::abort(paste0(
          "`PEP` contains finite values outside [0,1]. Examples: ",
          paste(ex, collapse = ", ")
        ))
      } else if (pep_mode == "sanitize") {
        n_bad <- length(bad_pep)
        pep[bad_pep] <- NA_real_
        rlang::warn(paste0(
          "PEP values outside [0,1] were set to NA (n=", n_bad, "). ",
          "Proceeding with q-values computed from Score + Reverse."
        ))
      }
    }
  }

  out <- tibble::tibble(
    id = as.character(x[["id"]]),
    run = as.character(x[["Raw file"]]),
    is_decoy = as.logical(is_decoy),
    pep = pep,
    score = score
  )

  if (anyNA(out$id) || any(out$id == "")) rlang::abort("Some `id` values are missing/empty.")
  if (anyNA(out$run) || any(out$run == "")) rlang::abort("Some `Raw file` values are missing/empty.")
  if (anyNA(out$is_decoy)) rlang::abort("`is_decoy` contains NA after parsing Reverse column.")
  if (!any(out$is_decoy)) rlang::abort("No decoys detected (Reverse never '+'); cannot compute q-values.")
  if (!any(!out$is_decoy)) rlang::abort("No targets detected (all Reverse '+'); cannot compute q-values.")

  out$q <- tdc_qvalues(score = out$score, is_decoy = out$is_decoy, add_decoy = add_decoy)

  q_source <- paste0("reconstructed_TDC_from_MaxQuant_Score(add_decoy=", add_decoy, ")")

  as_dfdr_tbl(
    out,
    unit = unit,
    scope = scope,
    q_source = q_source,
    q_max_export = NA_real_,
    provenance = provenance
  )
}

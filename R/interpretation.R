# Plain-language summaries of the p-value and PEP diagnostics, used by the
# HTML report. The text is plain ASCII (pi0, alpha, <=, >=, ~), so that it
# survives any session encoding; the report turns these tokens into symbols.

.u <- list(pi0 = "pi0", alpha = "alpha", le = "<=", ge = ">=",
           approx = "~", times = "x", minus = "-")

.fmt_n <- function(v) formatC(v, format = "d", big.mark = ",")
.fmt_p <- function(v) {
  if (!is.finite(v)) return("NA")
  if (v == 0) return("0")
  if (abs(v) < 1e-3) return(formatC(v, format = "e", digits = 2))
  formatC(v, format = "fg", digits = 3, flag = "#")
}
.fmt_pi0 <- function(v) if (is.finite(v)) formatC(v, format = "f", digits = 2) else "NA"
.fmt_range <- function(lo, hi) {
  if (!is.finite(lo) || !is.finite(hi)) return("NA")
  a <- .fmt_pi0(lo); b <- .fmt_pi0(hi)
  if (a == b) a else paste0(a, "-", b)
}

# Largest p-value threshold t such that the estimated FDR
# pi0 * m * t / #{p <= t} is at most alpha (Benjamini-Hochberg step-up when
# pi0 = 1; Storey-type when pi0 < 1). Returns NA when nothing is selected.
.stepup_threshold <- function(p, alpha, pi0 = 1) {
  p <- sort(p[is.finite(p)])
  m <- length(p)
  if (!m || !is.finite(pi0) || pi0 <= 0) return(NA_real_)
  ok <- which(p <= seq_len(m) * alpha / (m * min(pi0, 1)))
  if (!length(ok)) return(NA_real_)
  p[max(ok)]
}

.pi0_methods <- c("st.spline", "st.boot", "jiang", "histo", "langaas",
                  "pounds", "abh", "slim")

# pi0 by the eight cp4p estimators; falls back on Storey (lambda = 0.5) if
# cp4p fails (e.g. heavily tied p-values).
.pi0_all <- function(p, nbins = 20, pz = 0.05) {
  p <- p[is.finite(p) & p >= 0 & p <= 1]
  storey <- min(1, mean(p > 0.5) / 0.5)
  r <- tryCatch(cp4p::estim.pi0(p = p, pi0.method = "ALL", nbins = nbins, pz = pz),
                error = function(e) NULL)
  v <- if (!is.null(r)) suppressWarnings(as.numeric(r$pi0)) else numeric(0)
  if (length(v) == length(.pi0_methods) && any(is.finite(v))) {
    tibble::tibble(method = .pi0_methods, pi0 = pmin(1, pmax(0, v)), source = "cp4p")
  } else {
    tibble::tibble(method = "storey_0.5", pi0 = storey, source = "fallback")
  }
}

# Group the pi0 estimates: the estimators within `tol` of the median of all
# estimates converge; the consensus value is their median. With fewer than
# three converging estimators (out of three or more) there is no consensus and
# the median of all estimates is returned.
.pi0_consensus <- function(v, m, tol = 0.05) {
  ok <- is.finite(v)
  v <- v[ok]; m <- m[ok]
  if (!length(v)) return(list(value = NA_real_, n_in = 0L, members = character(0),
                              others = character(0), others_txt = "", lo = NA_real_,
                              hi = NA_real_, converged = FALSE))
  med <- stats::median(v)
  inn <- which(abs(v - med) <= tol)
  converged <- length(inn) >= min(3, length(v))
  idx <- if (converged) inn else seq_along(v)
  oth <- setdiff(seq_along(v), idx)
  list(value = stats::median(v[idx]), n_in = if (converged) length(inn) else 0L,
       members = m[idx], others = m[oth],
       others_txt = paste(sprintf("%s %s", m[oth], formatC(v[oth], format = "f", digits = 2)), collapse = ", "),
       lo = min(v[idx]), hi = max(v[idx]), converged = converged)
}

.p_kind <- function(x) {
  src <- (attr(x, "meta") %||% list())$p_source %||% NA_character_
  if (is.na(src)) return("supplied")
  if (grepl("-PEP", src, fixed = TRUE)) return("pep")
  if (grepl("score_to_pvalue(score", src, fixed = TRUE)) {
    if (grepl("decoy_ecdf", src, fixed = TRUE)) return("score_decoy_ecdf")
    return("score")
  }
  "supplied"
}

#' Summary of the p-value diagnostics: pi0 estimates and FDR thresholds
#'
#' Summarises what the p-value diagnostics of \code{\link{dfdr_run_all}} say
#' about a list: the proportion of true nulls \eqn{\pi_0} (incorrect
#' identifications) estimated on targets and decoys together, on targets and on
#' decoys, and the p-value threshold that gives an estimated FDR of
#' \code{alpha}, with the corresponding score (or PEP) threshold when the
#' p-values were derived from a score (or from \eqn{-}PEP).
#'
#' \eqn{\pi_0} is estimated with the eight estimators of
#' \code{cp4p::estim.pi0()} (Storey with spline or bootstrap, Jiang, histogram,
#' Langaas, Pounds, ABH, SLIM); if cp4p fails, Storey's estimator at
#' \eqn{\lambda = 0.5} is used. The estimates are then grouped: the estimators
#' within \code{tol} of the median of all estimates converge, and their median
#' is the consensus value; the others are listed separately. With fewer than
#' three converging estimators there is no consensus.
#'
#' The FDR threshold is the largest \eqn{t} such that
#' \eqn{\hat\pi_0 m t / \#\{p \le t\} \le \alpha} (the estimator of the
#' "Estimated FDR" plots), computed with the consensus \eqn{\hat\pi_0} (or the
#' estimator named in \code{pi0_method_fdr}) and with \eqn{\pi_0 = 1}
#' (Benjamini--Hochberg, the conservative choice).
#'
#' @param x A \code{dfdr_tbl} with a \code{p} column.
#' @param alpha FDR level (default 0.01).
#' @param nbins,pz Passed to \code{cp4p::estim.pi0()}.
#' @param pi0_method_fdr \eqn{\pi_0} used for the threshold: \code{"consensus"}
#'   (default) or the name of one estimator (e.g. \code{"pounds"}).
#' @param min_n Minimum number of p-values in a subset for \eqn{\pi_0} to be
#'   estimated (default 2000).
#' @param tol Maximum distance to the median of all estimates for an estimator
#'   to count as converging (default 0.05).
#'
#' @return A list with \code{pi0} (one row per subset and estimator),
#'   \code{summary} (one row per subset: count, the eight estimates, the
#'   consensus value, how many estimators converge and which, the others, and
#'   the \eqn{\pi_0} used for the threshold), \code{thresholds} (one row per
#'   subset and \eqn{\pi_0} choice: p-value threshold, number of accepted
#'   entries, targets and decoys, score or PEP threshold) and \code{info} (kind
#'   of p-values and their resolution).
#' @seealso \code{\link{dfdr_p_interpretation}}
#' @examples
#' set.seed(1)
#' n <- 1500
#' dec <- rep(c(TRUE, FALSE), c(500, 1000))
#' sc <- c(rnorm(500), rnorm(500), rnorm(500, 4))
#' df <- data.frame(id = paste0("u", seq_len(n)), is_decoy = dec, q = 1,
#'                  pep = NA_real_, run = "r", score = sc)
#' df$p <- score_to_pvalue(sc, method = "decoy_ecdf", is_decoy = dec)
#' x <- as_dfdr_tbl(df, unit = "psm", scope = "global", q_source = "toy")
#' attr(x, "meta")$p_source <- "score_to_pvalue(score, method='decoy_ecdf')"
#' if (requireNamespace("cp4p", quietly = TRUE)) {
#'   s <- dfdr_p_summary(x, alpha = 0.01, min_n = 400)
#'   s$summary
#'   s$thresholds
#' }
#' @export
dfdr_p_summary <- function(x, alpha = 0.01, nbins = 20, pz = 0.05,
                           pi0_method_fdr = "consensus", min_n = 2000, tol = 0.05) {
  if (!"p" %in% names(x) || !any(is.finite(x$p))) rlang::abort("No finite p-values in `x`.")
  kind <- .p_kind(x)
  ok <- is.finite(x$p) & x$p >= 0 & x$p <= 1
  dec <- x$is_decoy[ok]; p <- x$p[ok]
  sc <- if ("score" %in% names(x)) x$score[ok] else rep(NA_real_, sum(ok))
  pep <- if ("pep" %in% names(x)) x$pep[ok] else rep(NA_real_, sum(ok))
  nT <- sum(!dec); nD <- sum(dec)

  subsets <- list(all = rep(TRUE, length(p)), target = !dec, decoy = dec)
  pi0_tbl <- purrr::imap_dfr(subsets, function(sel, nm) {
    if (sum(sel) < min_n) return(tibble::tibble(subset = nm, n = sum(sel), method = NA_character_,
                                                pi0 = NA_real_, source = "too few p-values"))
    dplyr::mutate(.pi0_all(p[sel], nbins = nbins, pz = pz), subset = nm, n = sum(sel), .before = 1)
  })

  summ <- purrr::map_dfr(names(subsets), function(nm) {
    v <- pi0_tbl$pi0[pi0_tbl$subset == nm]
    m <- pi0_tbl$method[pi0_tbl$subset == nm]
    cons <- .pi0_consensus(v, m, tol = tol)
    used <- if (identical(pi0_method_fdr, "consensus")) cons$value else
      if (pi0_method_fdr %in% m) v[m == pi0_method_fdr][1] else cons$value
    wide <- stats::setNames(as.list(rep(NA_real_, length(.pi0_methods))), .pi0_methods)
    for (i in seq_along(m)) if (!is.na(m[i]) && m[i] %in% .pi0_methods) wide[[m[i]]] <- v[i]
    tibble::as_tibble(c(
      list(subset = nm, n = sum(subsets[[nm]])),
      wide,
      list(pi0_consensus = cons$value,
           n_converging = cons$n_in,
           n_estimators = sum(is.finite(v)),
           converging = paste(cons$members, collapse = ", "),
           others = cons$others_txt,
           pi0_used = used)))
  })

  thr <- purrr::map_dfr(c("target", "all"), function(nm) {
    sel <- subsets[[nm]]
    pi0_est <- summ$pi0_used[summ$subset == nm]
    purrr::map_dfr(list(c("estimate", pi0_est), c("1 (BH)", 1)), function(ch) {
      pi0v <- as.numeric(ch[2])
      t <- if (is.finite(pi0v)) .stepup_threshold(p[sel], alpha, pi0v) else NA_real_
      acc <- sel & is.finite(t) & p <= t
      accT <- acc & !dec
      s_thr <- if (any(accT) && kind %in% c("score", "score_decoy_ecdf") && any(is.finite(sc[accT])))
        min(sc[accT], na.rm = TRUE) else NA_real_
      pep_thr <- if (any(accT) && kind == "pep" && any(is.finite(pep[accT])))
        max(pep[accT], na.rm = TRUE) else NA_real_
      tibble::tibble(subset = nm, pi0_choice = ch[1], pi0 = pi0v, alpha = alpha,
                     p_threshold = t, n_accepted = sum(acc), n_accepted_target = sum(accT),
                     n_accepted_decoy = sum(acc & dec),
                     score_threshold = s_thr, pep_threshold = pep_thr)
    })
  })

  p_min_T <- if (any(!dec)) min(p[!dec]) else NA_real_
  list(pi0 = pi0_tbl, summary = summ, thresholds = thr,
              info = list(kind = kind, alpha = alpha, pi0_method_fdr = pi0_method_fdr,
                   two_pass = isTRUE((attr(x, "meta") %||% list())$two_pass),
                   q_max_export = (attr(x, "meta") %||% list())$q_max_export %||% NA_real_,
                   p_min_target = p_min_T,
                   resolution = if (kind %in% c("score_decoy_ecdf", "pep") && nD > 0) 1 / (nD + 1) else NA_real_,
                   p_source = (attr(x, "meta") %||% list())$p_source %||% NA_character_))
}

#' Plain-language interpretation of the p-value diagnostics
#'
#' Turns \code{\link{dfdr_p_summary}} into a few short statements for the
#' report: the \eqn{\pi_0} estimates (where the estimators converge, and which
#' differ), the p-value (and score or PEP) threshold giving an estimated FDR of
#' \code{alpha} on targets, compared with target--decoy counting, the pooled
#' targets + decoys check, and the caveats that apply. Cases handled: too few
#' p-values; decoy-calibrated pseudo-p-values (decoy \eqn{\pi_0 \approx 1} by
#' construction, so not a check) versus supplied p-values (decoy \eqn{\pi_0} is
#' then a real check); estimators that do not converge; no entry reaching the
#' FDR level; a threshold at the resolution limit of the pseudo-p-values;
#' two-pass exports.
#'
#' @param ps Output of \code{\link{dfdr_p_summary}}.
#' @param T_alpha Optional number of targets accepted by target--decoy
#'   counting (q \eqn{\le} alpha), for comparison.
#' @return A character vector of statements (plain ASCII; the report renders
#'   pi0, alpha, <= and >= as symbols).
#' @export
dfdr_p_interpretation <- function(ps, T_alpha = NULL) {
  u <- .u
  s <- ps$summary; th <- ps$thresholds; inf <- ps$info
  a_txt <- paste0(signif(inf$alpha * 100, 3), "%")
  row <- function(nm) s[s$subset == nm, , drop = FALSE]
  out <- character(0)
  caveats <- character(0)

  pi0_txt <- function(r, label) {
    if (!nrow(r) || !is.finite(r$pi0_consensus)) return(paste0(label, ": too few p-values"))
    ne <- r$n_estimators
    if (r$n_converging == 0) {
      v <- unlist(r[.pi0_methods]); v <- v[is.finite(v)]
      return(paste0(label, ": the estimators do not converge (", .fmt_range(min(v), max(v)), ")"))
    }
    if (ne == 1) return(paste0(label, " ", u$approx, " ", .fmt_pi0(r$pi0_consensus), " (Storey)"))
    v <- unlist(r[.pi0_methods]); mem <- strsplit(r$converging, ", ")[[1]]
    rng <- .fmt_range(min(v[mem]), max(v[mem]))
    paste0(label, " ", u$approx, " ", .fmt_pi0(r$pi0_consensus), " (",
           if (r$n_converging == ne) paste0("all ", ne, " estimators, ", rng)
           else paste0(r$n_converging, " of ", ne, " estimators converge, ", rng,
                       "; ", r$others), ")")
  }

  a <- row("all"); t <- row("target"); d <- row("decoy")
  if (!any(is.finite(s$pi0_consensus))) {
    return(paste0("Too few p-values to estimate ", u$pi0, " (fewer than the minimum per subset)."))
  }
  dec_txt <- if (nrow(d) && is.finite(d$pi0_consensus)) {
    if (inf$kind %in% c("score_decoy_ecdf", "pep")) paste0("decoys ", u$approx, " 1 by construction")
    else pi0_txt(d, "decoys")
  } else "decoys: too few p-values"
  out <- c(out, paste0("Proportion of incorrect identifications (", u$pi0, "): ",
                       pi0_txt(t, "targets"), "; ", pi0_txt(a, "targets + decoys"), "; ", dec_txt, "."))
  if (!inf$kind %in% c("score_decoy_ecdf", "pep") && nrow(d) && is.finite(d$pi0_consensus) &&
      abs(d$pi0_consensus - 1) > 0.05) {
    out <- c(out, paste0("Decoy p-values are not uniform (", u$pi0, " ", u$approx, " ", .fmt_pi0(d$pi0_consensus),
      "): the p-values are not calibrated or the decoys do not behave as incorrect matches. ",
      "After target-decoy competition the winner's p-value is not uniform even for an incorrect match; ",
      "use p-values computed before competition, or decoy-calibrated pseudo-p-values."))
  }

  # thresholds on targets
  tt <- th[th$subset == "target", , drop = FALSE]
  te <- tt[tt$pi0_choice == "estimate", , drop = FALSE]
  tb <- tt[tt$pi0_choice == "1 (BH)", , drop = FALSE]
  map_txt <- function(r) {
    if (is.finite(r$pep_threshold)) paste0(" (PEP ", u$le, " ", .fmt_p(r$pep_threshold), ")")
    else if (is.finite(r$score_threshold)) paste0(" (score ", u$ge, " ", formatC(r$score_threshold, format = "fg", digits = 4), ")")
    else ""
  }
  one <- function(r, lab) paste0("p ", u$le, " ", .fmt_p(r$p_threshold), map_txt(r), ", ",
                                 .fmt_n(r$n_accepted_target), " targets with ", lab)
  if (nrow(te) && !is.finite(te$p_threshold) && (!nrow(tb) || !is.finite(tb$p_threshold))) {
    out <- c(out, paste0("No target reaches an estimated FDR of ", a_txt, " (smallest target p-value ",
                         .fmt_p(inf$p_min_target), ")."))
  } else if (nrow(te)) {
    parts <- character(0)
    if (is.finite(te$p_threshold)) parts <- c(parts, one(te, paste0(u$pi0, " = ", .fmt_pi0(te$pi0))))
    if (nrow(tb) && is.finite(tb$p_threshold)) parts <- c(parts, one(tb, paste0(u$pi0, " = 1 (Benjamini-Hochberg)")))
    sent <- paste0("Estimated FDR of ", a_txt, " on targets: ", paste(parts, collapse = "; "), ".")
    if (!is.null(T_alpha) && is.finite(T_alpha)) {
      sent <- paste0(sent, " Target-decoy counting accepts ", .fmt_n(T_alpha), " targets at q ", u$le, " ", a_txt, ".")
    }
    out <- c(out, sent)
    if (is.finite(inf$resolution) && is.finite(te$p_threshold) && te$p_threshold <= inf$resolution * 1.0001) {
      caveats <- c(caveats, paste0("the threshold sits at the smallest attainable pseudo-p-value 1/(N_D + 1) = ",
                                   .fmt_p(inf$resolution), ", so the FDR level is not resolved more finely"))
    }
  }
  ta <- th[th$subset == "all" & th$pi0_choice == "estimate", , drop = FALSE]
  if (nrow(ta) && is.finite(ta$p_threshold)) {
    out <- c(out, paste0("Targets + decoys pooled (a check, not a list to report): p ", u$le, " ",
      .fmt_p(ta$p_threshold), " with ", u$pi0, " = ", .fmt_pi0(ta$pi0), ", ", .fmt_n(ta$n_accepted),
      " entries of which ", .fmt_n(ta$n_accepted_decoy), " decoys (",
      formatC(100 * ta$n_accepted_decoy / max(ta$n_accepted, 1), format = "f", digits = 2), "%)."))
  }

  if (inf$kind == "pep") caveats <- c(caveats, "pseudo-p-values from -PEP are coarse (tied PEPs; on lists collapsed across runs, each unit's smallest PEP)")
  if (isTRUE(inf$two_pass)) caveats <- c(caveats, paste0("two-pass export: the low-confidence region, from which ", u$pi0, " is estimated, is conditioned by the first pass"))
  if (length(caveats)) out <- c(out, paste0("Caveats: ", paste(caveats, collapse = "; "), "."))
  out
}

#' Summary of PEP concordance with the decoy-derived reference
#'
#' Compares, in PEP bins with mean target PEP \eqn{\le} \code{pep_max} and at
#' least \code{n_min} entries, the mean PEP of targets with the decoy-derived
#' reference \eqn{\hat e_b = D_b/T_b} (paper Eq. 1), and the sum of PEPs of the
#' targets accepted at \code{alpha} with the number of accepted decoys
#' \eqn{D_\alpha}, the target--decoy estimate of the number of false targets.
#'
#' @param x A \code{dfdr_tbl} with PEPs.
#' @param alpha Threshold for the accepted list (default 0.01).
#' @param binwidth,pep_max,n_min Binning, as in \code{\link{dfdr_ipe}}.
#' @param tol Signed differences smaller than \code{tol} (absolute) are read
#'   as agreement (default 0.01).
#' @return A list with \code{IPE}, \code{bins_used}, \code{ratio} (mean target PEP over mean reference), \code{signed} (weighted
#'   mean of mean target PEP minus reference), \code{share_under} (weight of
#'   bins where PEP is below the reference), \code{mean_pep}, \code{mean_ref},
#'   \code{sum_pep}, \code{D_alpha}, \code{T_alpha} and \code{verdict}
#'   (\code{"agree"}, \code{"under"}, \code{"over"}, \code{"mixed"},
#'   \code{"no_bins"}).
#' @seealso \code{\link{dfdr_pep_interpretation}}
#' @export
dfdr_pep_summary <- function(x, alpha = 0.01, binwidth = 0.05, pep_max = 0.5,
                             n_min = 200, tol = 0.01) {
  ipe <- dfdr_ipe(x, binwidth = binwidth, pep_max = pep_max, n_min = n_min)
  Ta <- sum(!x$is_decoy & x$q <= alpha); Da <- sum(x$is_decoy & x$q <= alpha)
  acc <- !x$is_decoy & x$q <= alpha & is.finite(x$pep)
  out <- list(IPE = ipe$IPE, bins_used = ipe$bins_used, signed = NA_real_, share_under = NA_real_,
              mean_pep = NA_real_, mean_ref = NA_real_, sum_pep = sum(x$pep[acc]),
              D_alpha = Da, T_alpha = Ta, alpha = alpha, verdict = "no_bins")
  b <- ipe$bins
  if (!nrow(b) || !any(b$used)) return(out)
  b <- b[b$used, ]
  w <- b$n / sum(b$n)
  d <- b$pep_mean_targets - b$e_b
  out$signed <- sum(w * d)
  out$share_under <- sum(w[d < -tol])
  share_over <- sum(w[d > tol])
  out$mean_pep <- sum(w * b$pep_mean_targets)
  out$mean_ref <- sum(w * b$e_b)
      out$ratio <- if (out$mean_ref > 0) out$mean_pep / out$mean_ref else NA_real_
  # Agreement needs both a small absolute gap (IPE) and a small relative one:
  # with small PEPs the IPE can stay below 0.05 while PEPs are several times
  # the reference.
  rel_ok <- is.finite(out$ratio) && out$ratio >= 0.8 && out$ratio <= 1.25
  out$verdict <- if (rel_ok && (!is.finite(out$IPE) || out$IPE <= 0.05)) "agree" else
    if (out$share_under >= 0.25 && share_over >= 0.25) "mixed" else
      if (out$signed < 0) "under" else "over"
  out
}

#' Plain-language interpretation of PEP concordance
#'
#' @param ps Output of \code{\link{dfdr_pep_summary}}.
#' @return A character vector of sentences.
#' @export
dfdr_pep_interpretation <- function(ps) {
  u <- .u
  a_txt <- paste0(signif(ps$alpha * 100, 3), "%")
  sum_txt <- if (ps$T_alpha > 0) paste0(" At ", u$alpha, " = ", a_txt, ", the PEPs of the ", .fmt_n(ps$T_alpha),
        " accepted targets sum to ", formatC(round(ps$sum_pep), format = "d", big.mark = ","), ", against ", .fmt_n(ps$D_alpha),
    " accepted decoys, the target", "-", "decoy estimate of the number of false targets") else ""
  ratio <- if (ps$D_alpha > 0 && ps$sum_pep > 0) ps$D_alpha / ps$sum_pep else NA_real_
    ref_txt <- paste0("(weighted mean target PEP ", .fmt_p(ps$mean_pep), " against a decoy-derived reference D_b/T_b of ",
                    .fmt_p(ps$mean_ref),
                    if (is.finite(ps$ratio %||% NA_real_)) paste0(", a ratio of ", formatC(ps$ratio, format = "f", digits = 2)) else "",
                    "; IPE = ", .fmt_p(ps$IPE),
                    if (is.finite(ps$IPE)) paste0(", ", if (ps$IPE > 0.05) "above" else "below", " the 0.05 alert threshold") else "",
                    ")")
  caveat <- " This is concordance with the decoys under equal chance, an internal check, not calibration against ground truth."
  switch(ps$verdict,
    no_bins = paste0("Too few entries per PEP bin (or no PEP below 0.5) to compare PEPs with the decoy-derived reference.", caveat),
            agree = paste0("The PEPs agree with the decoy-derived reference ", ref_txt,
                   if (is.finite(ps$signed) && ps$signed < -0.01) ", being slightly lower on average"
                   else if (is.finite(ps$signed) && ps$signed > 0.01) ", being slightly higher on average" else "",
                                                         ": read as probabilities, they match the error rate the decoys indicate.",
                   sum_txt, if (nzchar(sum_txt)) "." else "",
                   if (nzchar(sum_txt) && is.finite(ratio) && (ratio > 1.5 || ratio < 1 / 1.5))
                     paste0(" In the accepted list itself, where PEPs are smallest, the decoys indicate about ",
                            formatC(if (ratio > 1) ratio else 1 / ratio, format = "f", digits = 1), " times ",
                            if (ratio > 1) "more" else "fewer",
                            " false targets than the PEPs: the agreement holds on the PEP bins up to 0.5, less so at the top of the list.") else "",
                   caveat),
    under = paste0("The PEPs are lower than the decoy-derived reference ", ref_txt,
      ": they underestimate the probability that an identification is false, so reading them as probabilities overstates the confidence in the identifications (anti-conservative).",
      sum_txt, if (nzchar(sum_txt) && is.finite(ratio)) paste0(", about ", formatC(ratio, format = "f", digits = 1), " times more") else "", if (nzchar(sum_txt)) "." else "",
      " A common cause is PEPs computed for another unit (e.g. per run) and then aggregated, for instance by keeping each unit's smallest PEP across runs; use them as a ranking, or on the unit for which they were computed.", caveat),
    over = paste0("The PEPs are higher than the decoy-derived reference ", ref_txt,
      ": they overestimate the probability that an identification is false (conservative), so a PEP-based filter would discard correct identifications.",
      sum_txt, if (nzchar(sum_txt) && is.finite(ratio)) paste0(", about ", formatC(1 / ratio, format = "f", digits = 1), " times fewer") else "", if (nzchar(sum_txt)) "." else "", caveat),
    mixed = paste0("The PEPs are below the decoy-derived reference in some PEP bins and above it in others ", ref_txt,
      ": they rank identifications but their values are not reliable probabilities across the range.", sum_txt, if (nzchar(sum_txt)) "." else "", caveat))
}

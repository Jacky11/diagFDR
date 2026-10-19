#' Convert identification scores to p-values or pseudo-p-values
#'
#' Converts identification scores to:
#' \itemize{
#'   \item \strong{p-values} when the score has a known p-value definition (e.g. Mascot score),
#'   \item \strong{pseudo-p-values} when no such definition exists (e.g. MaxQuant/Andromeda scores).
#' }
#'
#' Pseudo-p-values are useful as inputs to plausibility checks (e.g. calibration/ECDF diagnostics),
#' but they are not guaranteed to be valid null p-values unless the method is explicitly null-based
#' (e.g. \code{method = "decoy_ecdf"} with a reliable decoy set).
#'
#' @param score Numeric vector. Identification scores. Higher-is-better is assumed for
#'   \code{method = "rank"} and \code{method = "decoy_ecdf"} (flip sign upstream if needed).
#' @param method Character. Conversion method:
#'   \describe{
#'     \item{\code{"mascot"}}{Treat \code{score} as Mascot score \eqn{S=-10\log_{10}(p)}; returns \eqn{p=10^{-S/10}}.}
#'     \item{\code{"neglog10p"}}{Generic \eqn{S=-10\log_{10}(p)}; same transform as \code{"mascot"}.}
#'     \item{\code{"evalue"}}{Treat \code{score} as an expectation/e-value-like quantity and return it as a p-like value.}
#'     \item{\code{"rank"}}{Rank-based pseudo-p-values \eqn{p_i=\mathrm{rank}(-s_i)/(n+1)}.}
#'     \item{\code{"decoy_ecdf"}}{Decoy-ECDF pseudo-p-values from the empirical right-tail of decoy scores;
#'       requires \code{is_decoy}.}
#'   }
#' @param is_decoy Optional logical vector (same length as \code{score}). Required for
#'   \code{method = "decoy_ecdf"}. TRUE=decoy, FALSE=target.
#' @param eps Numeric scalar > 0. Lower bound to avoid exact 0 (default \code{.Machine$double.xmin}).
#' @param clamp Logical. If \code{TRUE} (default), clamp output to \code{[0,1]}.
#' @param ties Character. Tie-handling for \code{method="rank"}; passed to \code{rank()}.
#'
#' @return
#' A numeric vector of p-values/pseudo-p-values of the same length as \code{score}.
#' Non-finite scores yield \code{NA_real_}. For \code{method="decoy_ecdf"}, returned
#' values estimate the decoy right-tail probability \eqn{P(S_{decoy}\ge s)} (with a
#' small +1 smoothing), which can be interpreted as a null-calibrated pseudo-p-value
#' when decoys provide a good null sample.
#'
#' @examples
#' set.seed(1)
#' score <- rnorm(10)
#'
#' # Mascot-like transformation: S = -10 log10(p)
#' S <- c(10, 20, 30)
#' score_to_pvalue(S, method = "mascot")
#'
#' # Rank-based pseudo-p-values (higher score => smaller p)
#' score_to_pvalue(score, method = "rank")
#'
#' # Decoy-ECDF pseudo-p-values (requires decoys)
#' n <- 200
#' score2 <- c(rnorm(n, mean = 0), rnorm(n, mean = 1))  # targets shifted higher
#' is_decoy <- c(rep(TRUE, n), rep(FALSE, n))
#' p2 <- score_to_pvalue(score2, method = "decoy_ecdf", is_decoy = is_decoy)
#' summary(p2)
#'
#' @export
score_to_pvalue <- function(score,
                            method = c("mascot", "neglog10p", "evalue", "rank", "decoy_ecdf"),
                            is_decoy = NULL,
                            eps = .Machine$double.xmin,
                            clamp = TRUE,
                            ties = c("average", "first", "random", "max", "min")) {
  method <- match.arg(method)
  ties <- match.arg(ties)

  if (!is.numeric(score)) rlang::abort("`score` must be numeric.")
  if (!is.numeric(eps) || length(eps) != 1 || !is.finite(eps) || eps <= 0) {
    rlang::abort("`eps` must be a single positive finite number.")
  }
  if (!is.logical(clamp) || length(clamp) != 1) rlang::abort("`clamp` must be TRUE/FALSE.")

  n <- length(score)
  ok <- is.finite(score)

  p <- rep(NA_real_, n)

  if (method %in% c("mascot", "neglog10p")) {
    # p = 10^(-S/10)
    p[ok] <- 10^(-score[ok] / 10)

  } else if (method == "evalue") {
    # p-like monotone mapping
    p[ok] <- score[ok]

  } else if (method == "rank") {
    # rank-based pseudo p-values, best score -> smallest p
    # rank(-score): largest score gets rank 1
    r <- rep(NA_real_, n)
    r[ok] <- rank(-score[ok], ties.method = ties)
    p[ok] <- r[ok] / (sum(ok) + 1)

  } else if (method == "decoy_ecdf") {
    if (is.null(is_decoy)) rlang::abort("`is_decoy` is required for method='decoy_ecdf'.")
    if (!is.logical(is_decoy) || length(is_decoy) != n) {
      rlang::abort("`is_decoy` must be a logical vector of the same length as `score`.")
    }
    if (anyNA(is_decoy)) rlang::abort("`is_decoy` contains NA; please fix upstream.")

    dec_scores <- score[ok & is_decoy]
    if (length(dec_scores) < 50) {
      rlang::abort("Too few finite decoy scores for method='decoy_ecdf' (need at least ~50).")
    }

    # Right-tail empirical survival function:
    # p(s) = P(S_decoy >= s) estimated from decoys
    # Use ranks among decoys to compute tail probability with +1 smoothing.
    dec_sorted <- sort(dec_scores, decreasing = FALSE) # ascending for findInterval
    # For each s, count decoys < s, then tail = (n_dec - count_lt + 1)/(n_dec + 1)
    count_lt <- findInterval(score[ok], dec_sorted, left.open = TRUE)  # number of decoys <=? approx
    n_dec <- length(dec_sorted)
    tail <- (n_dec - count_lt + 1) / (n_dec + 1)
    p[ok] <- tail
  }

  # sanitize numeric issues
  p[is.finite(p)] <- pmax(p[is.finite(p)], eps)
  if (isTRUE(clamp)) p[is.finite(p)] <- pmin(p[is.finite(p)], 1)

  p
}



#' Add pseudo p-values to \code{dfdr_tbl} objects that have scores but no p-values
#'
#' Iterates over a named list of \code{dfdr_tbl} objects and, where a finite \code{p}
#' column is missing but \code{score} and \code{is_decoy} are available, adds a new
#' \code{p} column computed by \code{\link{score_to_pvalue}}. Updates the metadata
#' attribute \code{p_source}.
#'
#' @param xs Named list of \code{dfdr_tbl} objects.
#' @param method Character. Method passed to \code{\link{score_to_pvalue}} (default \code{"decoy_ecdf"}).
#' @param score_source Character. \code{"score"} (default) uses the \code{score}
#'   column only. \code{"pep"} uses \code{score} when it is available and
#'   otherwise \code{-pep} as a monotone surrogate score (see
#'   \code{\link{dfdr_run_all}} for the caveats).
#'
#' @return
#' A list of \code{dfdr_tbl} objects of the same length/order as \code{xs}, with a
#' \code{p} column added where applicable. When no p-value can be computed, the
#' list is returned unchanged, its metadata field \code{p_source} records why
#' and a warning names the list.
#'
#' @keywords internal
add_pseudo_pvalues_to_list <- function(xs, method = "decoy_ecdf",
                                       score_source = c("score", "pep")) {
  score_source <- match.arg(score_source)
  out <- purrr::imap(xs, function(x, nm) {
    meta <- attr(x, "meta") %||% list()

    # Keep p-values that are already there
    if ("p" %in% names(x) && any(is.finite(x$p))) {
      return(x)
    }

    has_score <- "score" %in% names(x) && any(is.finite(x$score))
    has_pep <- "pep" %in% names(x) && any(is.finite(x$pep))

    if (has_score) {
      s <- x$score
      src <- sprintf("score_to_pvalue(score, method='%s')", method)
    } else if (score_source == "pep" && has_pep) {
      s <- -x$pep
      src <- sprintf("score_to_pvalue(-PEP, method='%s'; PEP used as surrogate score)", method)
    } else {
      why <- if (score_source == "pep") "no score and no PEP column" else "no score column"
      meta$p_source <- sprintf("not computed (%s)", why)
      attr(x, "meta") <- meta
      rlang::warn(sprintf(paste0(
        "List '%s': pseudo-p-values not computed (%s). ",
        "Pass `score_col` to the import function if the export has a discriminant score ",
        "(DIA-NN 2.x reports have no CScore column), or use ",
        "`pseudo_pvalue_score = \"pep\"` in dfdr_run_all() to use -PEP as a surrogate score."),
        nm, why))
      return(x)
    }

    p <- tryCatch(
      score_to_pvalue(score = s, method = method, is_decoy = x$is_decoy),
      error = function(e) {
        rlang::warn(sprintf("List '%s': pseudo-p-values not computed (%s).", nm, conditionMessage(e)))
        NULL
      })
    if (is.null(p)) {
      meta$p_source <- "not computed (score_to_pvalue() failed)"
      attr(x, "meta") <- meta
      return(x)
    }

    x$p <- p
    meta$p_source <- src
    attr(x, "meta") <- meta
    x
  })
  names(out) <- names(xs)
  out
}

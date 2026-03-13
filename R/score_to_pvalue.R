#' Convert identification scores to p-values or pseudo-p-values
#'
#' Converts identification scores to:
#' \itemize{
#'   \item \strong{p-values} when the score has a known p-value definition (e.g. Mascot score),
#'   \item \strong{pseudo-p-values} when no such definition exists (e.g. MaxQuant/Andromeda scores).
#' }
#'
#' Pseudo-p-values are useful as inputs to calibration checks (e.g. CP4P) and for building
#' common scales across tools, but they are not guaranteed to be valid null p-values unless
#' the method is explicitly null-based (e.g. \code{method="decoy_ecdf"} with reliable decoys).
#'
#' @param score Numeric vector. Identification score(s). Higher-is-better is assumed for
#'   \code{method="rank"} and \code{method="decoy_ecdf"} (you can flip sign upstream if needed).
#' @param method Character scalar. One of:
#'   \describe{
#'     \item{\code{"mascot"}}{Treat \code{score} as Mascot score: \eqn{S = -10\log_{10}(p)}.
#'       Returns \eqn{p = 10^{-S/10}}.}
#'     \item{\code{"neglog10p"}}{Generic \eqn{S = -10\log_{10}(p)}. Same transform as \code{"mascot"}.}
#'     \item{\code{"evalue"}}{Treat \code{score} as an E-value/expectation-like quantity.
#'       Returns \eqn{p = \min(1, \max(score, eps))}. This is a pseudo-p-value unless you have
#'       additional guarantees.}
#'     \item{\code{"rank"}}{Pseudo-p-values from score ranks:
#'       \eqn{p_i = \frac{\mathrm{rank}(-s_i)}{n+1}} (best score gets smallest p).
#'       This is monotone and in (0,1), but not a null p-value.}
#'     \item{\code{"decoy_ecdf"}}{Pseudo-p-values from empirical decoy null tail.
#'       Requires \code{is_decoy}. For each score \eqn{s}, returns
#'       \eqn{p(s) = P(S_{decoy} \ge s)} estimated from the decoy score ECDF (right-tail).
#'       This is the most defensible pseudo-p-value for arbitrary scores when decoys are reliable.}
#'   }
#' @param is_decoy Optional logical vector (same length as \code{score}). Required for
#'   \code{method="decoy_ecdf"}. TRUE=decoy, FALSE=target.
#' @param eps Numeric scalar > 0. Lower bound to avoid exact 0. Default \code{.Machine$double.xmin}.
#' @param clamp Logical. If TRUE (default), clamp output to `[0,1]`.
#' @param ties Character. How to handle ties for \code{method="rank"}. Passed to \code{rank()}.
#'   Default \code{"average"}.
#'
#' @return Numeric vector of p-values/pseudo-p-values (same length as \code{score}).
#'   Missing/non-finite scores yield \code{NA_real_}.
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



#' Add pseudo p-values to dfdr_tbl objects that have scores but no p-values
#'
#' @param xs Named list of dfdr_tbl objects
#' @param method Method passed to score_to_pvalue (default "decoy_ecdf")
#' @return Modified list with p column added where applicable
#' @keywords internal
add_pseudo_pvalues_to_list <- function(xs, method = "decoy_ecdf") {
  purrr::map(xs, function(x) {
    # Skip if already has p-values
    if ("p" %in% names(x) && any(is.finite(x$p))) {
      return(x)
    }

    # Skip if no score column
    if (!"score" %in% names(x) || all(is.na(x$score))) {
      return(x)
    }

    # Compute pseudo p-values
    x$p <- score_to_pvalue(
      score = x$score,
      method = method,
      is_decoy = x$is_decoy
    )

    # Update metadata
    meta <- attr(x, "meta") %||% list()
    meta$p_source <- sprintf("score_to_pvalue(method='%s')", method)
    attr(x, "meta") <- meta

    x
  })
}

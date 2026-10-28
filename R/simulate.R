## =====================================================================
##  Ground-truth simulation of target--decoy competition
##  (generative model of Supplementary Note S6 of the diagFDR paper)
## =====================================================================

#' Simulate a competed target/decoy list with known ground truth
#'
#' Each of \code{n} units carries a target score \eqn{X_i} and a decoy score
#' \eqn{Y_i}. A fraction \code{pi0} of units are null, with both scores drawn
#' independently from the same standardised null \eqn{F_0}; the remainder
#' carry a correct match, with the target score shifted by \code{delta}.
#' Under competition the unit enters the list with score \eqn{\max(X_i,Y_i)}
#' and label target if \eqn{X_i>Y_i}, decoy otherwise, so equal chance holds
#' by construction for null units. Q-values are computed by target--decoy
#' counting (\eqn{D/T}, running minimum), and the true PEP of each
#' competed target is returned for Gaussian nulls. Because the null status of
#' every unit is recorded, the realised FDP of any accepted list is known
#' exactly (\code{\link{dfdr_realised_fdp}}).
#'
#' Decoy misspecification can be simulated with \code{decoy_shift}
#' (decoys too easy to beat, scenario S1 of the paper) and
#' \code{decoy_scale} (decoys under-dispersed, S2); \code{compete = FALSE}
#' lists targets and decoys separately (S3); \code{rho} correlates scores
#' within blocks of \code{block} units (S5).
#'
#' @param n Number of units.
#' @param pi0 Proportion of null units.
#' @param delta Location shift of correct target scores.
#' @param null Null family: \code{"gaussian"}, \code{"gumbel"} or
#'   \code{"student"} (\eqn{t_4}), each standardised to mean 0 and variance 1.
#' @param decoy_shift Added to every decoy score.
#' @param decoy_scale Multiplies every decoy score.
#' @param compete Pair each target with its decoy (\code{TRUE}) or list
#'   both (\code{FALSE}).
#' @param pi0_lib With \code{compete = FALSE}, the null fraction among
#'   listed targets (an enriched library).
#' @param rho,block Within-block correlation and block size (Gaussian
#'   nulls only; \code{rho = 0} for independence).
#' @param add_decoy Added to the decoy count in the q-value (0 or 1).
#' @param seed Optional seed.
#'
#' @return A \code{dfdr_tbl} with columns \code{id}, \code{is_decoy},
#'   \code{q}, \code{pep} (true PEP, Gaussian competed lists only; otherwise
#'   \code{NA}), \code{run} (\code{NA}) and \code{score}. The attribute
#'   \code{"truth"} is a data.frame with \code{is_null} (the unit is
#'   incorrect) aligned with the rows, and the attribute \code{"params"}
#'   records the design.
#' @examples
#' sim <- dfdr_simulate_tdc(n = 20000, seed = 1)
#' dfdr_headline(sim, alpha = 0.01)
#' dfdr_realised_fdp(sim, alpha = c(0.01, 0.05))
#' @export
dfdr_simulate_tdc <- function(n = 2e5, pi0 = 0.667, delta = 4.65,
                              null = c("gaussian", "gumbel", "student"),
                              decoy_shift = 0, decoy_scale = 1,
                              compete = TRUE, pi0_lib = pi0,
                              rho = 0, block = 10L, add_decoy = 0L, seed = NULL) {
  null <- match.arg(null)
  if (!is.null(seed)) set.seed(seed)
  rnull <- switch(null,
    gaussian = function(m) stats::rnorm(m),
    gumbel   = function(m) { g <- -log(-log(stats::runif(m))); (g - 0.5772157) / (pi / sqrt(6)) },
    student  = function(m) stats::rt(m, df = 4) / sqrt(2))

  if (compete) {
    is_null <- stats::runif(n) < pi0
    if (rho > 0) {
      if (null != "gaussian") rlang::abort("rho > 0 is implemented for Gaussian nulls only.")
      g  <- rep(seq_len(ceiling(n / block)), each = block)[seq_len(n)]
      b  <- stats::rnorm(max(g)); cc <- stats::rnorm(max(g))
      X <- sqrt(rho) * b[g] + sqrt(1 - rho) * stats::rnorm(n)
      Y <- sqrt(rho) * cc[g] + sqrt(1 - rho) * stats::rnorm(n)
    } else {
      X <- rnull(n); Y <- rnull(n)
    }
    X[!is_null] <- X[!is_null] + delta
    Y <- Y * decoy_scale + decoy_shift
    is_decoy <- Y > X
    score <- pmax(X, Y)
    truth <- data.frame(is_null = is_null)
    ## true PEP of a competed target under Gaussian independent nulls:
    ## P(null | target wins with score s) with f_T(s) = (1-pi0) phi(s-delta) Phi(s) + pi0 phi(s) Phi(s)
    pep <- if (null == "gaussian" && rho == 0 && decoy_shift == 0 && decoy_scale == 1) {
      s <- score
      num <- pi0 * stats::dnorm(s)
      den <- num + (1 - pi0) * stats::dnorm(s - delta)
      p <- num / den
      p[is_decoy] <- NA_real_
      p
    } else NA_real_
  } else {
    ## separate listing: n targets with null fraction pi0_lib, n decoys
    is_null_t <- stats::runif(n) < pi0_lib
    X <- rnull(n); X[!is_null_t] <- X[!is_null_t] + delta
    Y <- rnull(n) * decoy_scale + decoy_shift
    score <- c(X, Y)
    is_decoy <- c(rep(FALSE, n), rep(TRUE, n))
    truth <- data.frame(is_null = c(is_null_t, rep(TRUE, n)))
    pep <- NA_real_
  }

  q <- tdc_qvalues(score, is_decoy, add_decoy = add_decoy)
  out <- data.frame(id = sprintf("u%07d", seq_along(score)), is_decoy = is_decoy,
                    q = q, pep = pep, run = NA_character_, score = score,
                    stringsAsFactors = FALSE)
  tbl <- as_dfdr_tbl(out, unit = "simulated_unit", scope = "global",
                     q_source = "tdc_counting", q_max_export = 1,
                     provenance = list(tool = "dfdr_simulate_tdc", null = null,
                                       pi0 = pi0, delta = delta, compete = compete),
                     two_pass = FALSE)
  attr(tbl, "truth") <- truth
  attr(tbl, "params") <- list(n = n, pi0 = pi0, delta = delta, null = null,
                              decoy_shift = decoy_shift, decoy_scale = decoy_scale,
                              compete = compete, pi0_lib = pi0_lib, rho = rho,
                              block = block, add_decoy = add_decoy)
  tbl
}

#' Realised FDP of a simulated list
#'
#' @param sim A list from \code{\link{dfdr_simulate_tdc}}.
#' @param alpha Threshold(s).
#' @return A data.frame with \code{alpha}, \code{T_alpha}, \code{D_alpha},
#'   \code{V_alpha} (accepted incorrect targets), \code{decoy_ratio} and the
#'   realised \code{FDP}.
#' @examples
#' sim <- dfdr_simulate_tdc(n = 20000, seed = 3)
#' dfdr_realised_fdp(sim, alpha = 0.01)
#' @export
dfdr_realised_fdp <- function(sim, alpha = 0.01) {
  truth <- attr(sim, "truth")
  if (is.null(truth)) rlang::abort("`sim` has no \"truth\" attribute; use dfdr_simulate_tdc().")
  do.call(rbind, lapply(alpha, function(a) {
    acc <- sim$q <= a
    Ta <- sum(acc & !sim$is_decoy); Da <- sum(acc & sim$is_decoy)
    Va <- sum(acc & !sim$is_decoy & truth$is_null)
    data.frame(alpha = a, T_alpha = Ta, D_alpha = Da, V_alpha = Va,
               decoy_ratio = if (Ta > 0) Da / Ta else NA_real_,
               FDP = if (Ta > 0) Va / Ta else NA_real_)
  }))
}

#' Monte Carlo check of the FDP-CV benchmark
#'
#' Repeats \code{\link{dfdr_simulate_tdc}} \code{B} times and compares the
#' coefficient of variation of the realised FDP across replicates with the
#' benchmark \eqn{\sqrt2\,D_\alpha^{-1/2}} of Proposition 1 of the diagFDR
#' paper, returning the ratio \eqn{\kappa = \mathrm{CV(FDP)} /
#' \overline{D_\alpha^{-1/2}}} (expected \eqn{\approx 1.41} under
#' independent Gaussian nulls).
#'
#' @param B Number of replicates.
#' @param alpha Threshold(s).
#' @param ... Passed to \code{\link{dfdr_simulate_tdc}} (do not pass
#'   \code{seed}).
#' @return A data.frame with, per \code{alpha}: mean FDP, \code{CV_FDP},
#'   mean \code{D_alpha}, mean support index, \code{kappa}, the share of
#'   replicates with FDP above \code{alpha}, \code{Var_V} and
#'   \code{2 E[D] + Var(D)}.
#' @examples
#' dfdr_check_cv_benchmark(B = 20, n = 20000, alpha = 0.01)
#' @export
dfdr_check_cv_benchmark <- function(B = 200, alpha = 0.01, ...) {
  res <- lapply(seq_len(B), function(b) dfdr_realised_fdp(dfdr_simulate_tdc(...), alpha = alpha))
  res <- do.call(rbind, res)
  do.call(rbind, lapply(alpha, function(a) {
    r <- res[res$alpha == a, ]
    cv <- stats::sd(r$FDP) / mean(r$FDP)
    si <- mean(r$D_alpha^(-0.5), na.rm = TRUE)
    data.frame(alpha = a, B = B, mean_FDP = mean(r$FDP), CV_FDP = cv,
               mean_D = mean(r$D_alpha), mean_support_index = si,
               kappa = cv / si, benchmark = sqrt(2) * si,
               share_FDP_above_alpha = mean(r$FDP > a),
               Var_V = stats::var(r$V_alpha),
               two_ED_plus_VarD = 2 * mean(r$D_alpha) + stats::var(r$D_alpha))
  }))
}

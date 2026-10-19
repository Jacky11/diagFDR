#' Evidential basis of each review alert
#'
#' Returns, for each flag produced by \code{\link{flag_headline}}, the
#' diagnostic it applies to, its alert threshold and the evidential basis of
#' that threshold, so that alerts are not read as equally calibrated:
#' \describe{
#'   \item{derived check}{follows from a stated result (e.g. the recomputed
#'     decoy ratio reproduces the definition of a counted q-value);}
#'   \item{simulation-informed alert}{threshold calibrated against the
#'     ground-truth simulations of the diagFDR paper;}
#'   \item{heuristic alert}{practical guideline taken from the paper's case
#'     studies, to be re-derived locally;}
#'   \item{descriptive}{reported without threshold, as a sensitivity measure.}
#' }
#' Entrapment-based estimates, when supplied, are external evidence and carry
#' no package alert. The decoy-ratio alert levels (1.2 and 1.5 times alpha)
#' are review margins around a derived identity: for q-values counted on the
#' reported list, \eqn{D_\alpha/T_\alpha \le \alpha} by construction.
#'
#' @return A \link[tibble:tibble]{tibble} with columns \code{flag},
#'   \code{diagnostic}, \code{threshold} and \code{basis}.
#' @examples
#' dfdr_alert_basis()
#' @export
dfdr_alert_basis <- function() {
  tibble::tibble(
    flag = c("flag_FDR", "flag_Dalpha", "flag_support", "flag_Dwin",
             "flag_IPE", "flag_equalchance", "jaccard / elasticity"),
    diagnostic = c("Recomputed decoy ratio D/T on the reported list",
                   "Accepted decoy count D_alpha",
                   "Decoy-support index D_alpha^-1/2",
                   "Local boundary support D_alpha,win",
                   "Decoy-reference PEP concordance (IPE)",
                   "Equal chance under paired target-decoy competition",
                   "List overlap and threshold elasticity"),
    threshold = c("> 1.5 alpha (caution > 1.2 alpha); expected <= alpha by construction", "< 10 (caution < 50)",
                  "> 0.30 (caution > 0.15)", "< 5 (caution < 20)", "> 0.05",
                  "|pi_D - 0.5| > 0.15 (caution > 0.10); n/e on two-pass exports",
                  "none"),
    basis = c("derived check", "simulation-informed alert",
              "simulation-informed alert", "heuristic alert", "heuristic alert",
              "simulation-informed alert", "descriptive")
  )
}

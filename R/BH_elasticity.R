#' BH list elasticity (Jaccard) across alpha values
#'
#' Computes the stability (elasticity) of Benjamini--Hochberg (BH) discovery sets
#' under a multiplicative perturbation of the nominal FDR level. For each
#' \code{alpha} in \code{alphas}, the function compares the BH discovery set at
#' \code{alpha} to the discovery set at \code{(1+eps)*alpha} using the Jaccard
#' index.
#'
#' @param x A \code{dfdr_tbl} (or data frame) with columns \code{id} and \code{p}.
#'   Only finite \code{p}-values are used.
#' @param alphas Numeric vector of BH levels in \eqn{(0,1)}.
#' @param eps Numeric scalar \eqn{\epsilon > 0}. Relative perturbation (default
#'   0.2), i.e. compares \code{alpha} vs \code{(1+eps)*alpha}.
#'
#' @return
#' A \link[tibble:tibble]{tibble} with one row per element of \code{alphas}. Columns include:
#' \describe{
#'   \item{alpha}{The nominal BH level.}
#'   \item{alpha_perturbed}{The perturbed level \code{min((1+eps)*alpha, 1)}.}
#'   \item{t_alpha}{BH rejection threshold at \code{alpha} (NA if no rejections).}
#'   \item{t_alpha_perturbed}{BH rejection threshold at \code{alpha_perturbed} (NA if none).}
#'   \item{R_alpha}{Number of discoveries at \code{alpha}.}
#'   \item{R_alpha_perturbed}{Number of discoveries at \code{alpha_perturbed}.}
#'   \item{jaccard}{Jaccard similarity between the two discovery ID sets.}
#' }
#' This output is intended to quantify how sensitive BH discoveries are to small
#' changes in the chosen FDR level.
#'
#' @examples
#' library(tibble)
#'
#' set.seed(1)
#' n <- 5000
#' x <- tibble(
#'   id = as.character(seq_len(n)),
#'   # mixture: mostly null p-values + a small enriched component
#'   p = c(stats::runif(4500), stats::rbeta(500, 0.3, 1))
#' )
#'
#' dfdr_bh_elasticity(x, alphas = c(1e-3, 5e-3, 1e-2, 2e-2), eps = 0.2)
#'
#' @export
dfdr_bh_elasticity <- function(x,
                               alphas = c(1e-04, 5e-04, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05),
                               eps = 0.2) {
  x <- .safe_df(x)
  .check_has_cols(x, c("id", "p"))

  if (!is.numeric(alphas) || any(!is.finite(alphas)) || any(alphas <= 0) || any(alphas >= 1)) {
    rlang::abort("`alphas` must be finite values in (0,1).")
  }
  if (!is.numeric(eps) || length(eps) != 1 || !is.finite(eps) || eps <= 0) {
    rlang::abort("`eps` must be a single positive number.")
  }

  ok <- is.finite(x$p)
  p <- x$p[ok]
  ids <- x$id[ok]

  out <- vector("list", length(alphas))

  for (i in seq_along(alphas)) {
    a <- alphas[i]
    a2 <- min((1 + eps) * a, 1 - 1e-12)

    t1 <- .bh_threshold(p, a)
    t2 <- .bh_threshold(p, a2)

    A1 <- if (is.finite(t1)) ids[p <= t1] else character()
    A2 <- if (is.finite(t2)) ids[p <= t2] else character()

    out[[i]] <- tibble::tibble(
      alpha = a,
      alpha_perturbed = a2,
      t_alpha = if (is.finite(t1)) t1 else NA_real_,
      t_alpha_perturbed = if (is.finite(t2)) t2 else NA_real_,
      R_alpha = as.integer(length(A1)),
      R_alpha_perturbed = as.integer(length(A2)),
      jaccard = .jaccard(A1, A2)
    )
  }

  dplyr::bind_rows(out)
}


#' Plot BH elasticity (Jaccard) vs alpha
#'
#' Convenience plotting function for the output of \code{\link{dfdr_bh_elasticity}}.
#' The x-axis is \code{log10(alpha)} and the y-axis is the Jaccard overlap between
#' BH discovery sets at \code{alpha} and \code{(1+eps)*alpha}.
#'
#' @param el_tbl A tibble as returned by \code{\link{dfdr_bh_elasticity}}.
#' @param xlab Character. X-axis label.
#' @param title Character. Plot title.
#'
#' @return
#' A \link[ggplot2:ggplot]{ggplot} object.
#'
#' @examples
#' library(tibble)
#'
#' set.seed(1)
#' n <- 2000
#' x <- tibble(
#'   id = as.character(seq_len(n)),
#'   p = c(stats::runif(1800), stats::rbeta(200, 0.3, 1))
#' )
#'
#' el <- dfdr_bh_elasticity(x, alphas = c(1e-3, 5e-3, 1e-2, 2e-2), eps = 0.2)
#' p <- dfdr_plot_bh_elasticity(el)
#' p
#'
#' @export
dfdr_plot_bh_elasticity <- function(el_tbl,
                                    xlab = "alpha (log10)",
                                    title = "BH elasticity: Jaccard(alpha, (1+eps)alpha)") {
  .check_has_cols(el_tbl, c("alpha", "jaccard"))
  ggplot2::ggplot(el_tbl, ggplot2::aes(x = log10(alpha), y = jaccard)) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1.2) +
    ggplot2::labs(x = xlab, y = "Jaccard", title = title) +
    ggplot2::ylim(0, 1)
}


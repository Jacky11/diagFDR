#' BH list elasticity (Jaccard) across alpha values
#'
#' Computes Jaccard overlap of BH discovery sets between alpha and (1+eps)*alpha.
#'
#' @param x A dfdr_tbl with columns \code{id} and \code{p}.
#' @param alphas Numeric vector of BH levels in (0,1).
#' @param eps Relative perturbation (default 0.2).
#'
#' @return Tibble with columns \code{alpha}, \code{alpha_perturbed}, \code{jaccard},
#'   \code{R_alpha}, \code{R_alpha_perturbed}, and thresholds.
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
#' @param el_tbl Output of \code{\link{dfdr_bh_elasticity}}.
#' @param xlab X-axis label.
#' @param title Plot title.
#' @return A ggplot object.
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

#' Storey \eqn{\pi_0(\lambda)} tail-uniformity diagnostic
#'
#' Estimates \eqn{\pi_0(\lambda) = \#\{p>\lambda\} / ((1-\lambda)m)} over a grid
#' of \code{lambdas}. This is a plausibility diagnostic for (pseudo-)p-values:
#' under a well-calibrated null, the upper tail should be approximately uniform,
#' leading to stable \eqn{\pi_0(\lambda)} curves.
#'
#' @param x A \code{dfdr_tbl} (or data.frame) containing column \code{p}.
#' @param lambdas Numeric vector in \eqn{[0,1)}. Typical range: 0.5--0.95.
#' @param stratify Optional character vector of column names used to stratify the
#'   diagnostic (e.g. \code{c("run")}).
#' @param min_n Integer. Minimum number of finite p-values required per stratum.
#' @param clamp Logical. If \code{TRUE} (default), clamp \code{pi0_hat} into
#'   \code{[0,1]} for reporting.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{pi0}{A \link[tibble:tibble]{tibble} with one row per \code{lambda} per
#'   stratum and columns \code{stratum}, \code{lambda}, \code{pi0_hat}, and \code{n}
#'   (the number of finite p-values in the stratum).}
#'   \item{summary}{A tibble with one row per stratum, including \code{n} and
#'   summary statistics of \code{pi0_hat} across \code{lambdas}
#'   (\code{median_pi0}, \code{sd_pi0}, \code{iqr_pi0}). Strata with \code{n < min_n}
#'   are reported with \code{NA} metrics and a \code{note}.}
#' }
#'
#' @examples
#' library(tibble)
#'
#' set.seed(1)
#' n <- 6000
#' df <- tibble(
#'   run = sample(c("run1", "run2"), n, replace = TRUE),
#'   # mostly uniform p-values + a small enriched component
#'   p = c(stats::runif(5400), stats::rbeta(600, 0.3, 1))
#' )
#'
#' out <- dfdr_pi0_storey(df, stratify = "run", lambdas = seq(0.5, 0.9, by = 0.1), min_n = 500)
#' head(out$pi0)
#' out$summary
#'
#' @export
dfdr_pi0_storey <- function(x,
                            lambdas = seq(0.5, 0.95, by = 0.05),
                            stratify = NULL,
                            min_n = 2000,
                            clamp = TRUE) {
  x <- .safe_df(x)
  .check_has_cols(x, c("p"))

  if (!is.numeric(lambdas) || any(!is.finite(lambdas)) || any(lambdas < 0) || any(lambdas >= 1)) {
    rlang::abort("`lambdas` must be numeric values in [0,1).")
  }
  lambdas <- sort(unique(lambdas))

  strata_df <- .get_strata_df(x, stratify)
  if (is.null(strata_df)) {
    x$stratum <- "all"
  } else {
    x$stratum <- .make_stratum_id(strata_df)
  }

  out_pi0 <- list()
  out_sum <- list()

  for (s in unique(x$stratum)) {
    xs <- x[x$stratum == s, , drop = FALSE]
    p <- xs$p
    p <- p[is.finite(p)]
    m <- length(p)

    if (m < min_n) {
      out_sum[[s]] <- tibble::tibble(
        stratum = s, n = m,
        median_pi0 = NA_real_,
        sd_pi0 = NA_real_,
        iqr_pi0 = NA_real_,
        note = paste0("Skipped: n < ", min_n)
      )
      next
    }

    pi0_hat <- vapply(lambdas, function(lam) {
      sum(p > lam) / ((1 - lam) * m)
    }, numeric(1))

    if (isTRUE(clamp)) pi0_hat <- pmin(pmax(pi0_hat, 0), 1)

    pi0_tbl <- tibble::tibble(
      stratum = s,
      lambda = lambdas,
      pi0_hat = pi0_hat,
      n = m
    )
    out_pi0[[s]] <- pi0_tbl

    out_sum[[s]] <- tibble::tibble(
      stratum = s,
      n = m,
      median_pi0 = stats::median(pi0_hat, na.rm = TRUE),
      sd_pi0 = stats::sd(pi0_hat, na.rm = TRUE),
      iqr_pi0 = stats::IQR(pi0_hat, na.rm = TRUE),
      note = NA_character_
    )
  }

  pi0_all <- if (length(out_pi0)) dplyr::bind_rows(out_pi0) else tibble::tibble()
  sum_all <- dplyr::bind_rows(out_sum)

  list(pi0 = pi0_all, summary = sum_all)
}

#' Plot Storey \eqn{\pi_0(\lambda)} curve
#'
#' Plots \eqn{\pi_0(\lambda)} estimates returned by \code{\link{dfdr_pi0_storey}}.
#' Values closer to 1 suggest many nulls (few true effects), while values far below
#' 1 (especially if unstable across \code{lambda}) may indicate deviations from
#' tail-uniformity.
#'
#' @param pi0_tbl A tibble, typically \code{out$pi0} from \code{\link{dfdr_pi0_storey}}.
#'   Must contain columns \code{stratum}, \code{lambda}, and \code{pi0_hat}.
#' @param title Character. Plot title.
#'
#' @return
#' A \link[ggplot2:ggplot]{ggplot} object.
#'
#' @examples
#' library(tibble)
#'
#' set.seed(1)
#' df <- tibble(
#'   stratum = rep("all", 5),
#'   lambda = seq(0.5, 0.9, by = 0.1),
#'   pi0_hat = c(0.95, 0.96, 0.97, 0.98, 0.99),
#'   n = 1000
#' )
#' p <- dfdr_plot_pi0(df)
#' p
#'
#' @export
dfdr_plot_pi0 <- function(pi0_tbl, title = "Storey pi0(lambda) diagnostic") {
  .check_has_cols(pi0_tbl, c("stratum", "lambda", "pi0_hat"))
  ggplot2::ggplot(pi0_tbl, ggplot2::aes(x = .data$lambda, y = .data$pi0_hat)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::facet_wrap(~ stratum, scales = "free_y") +
    ggplot2::labs(x = expression(lambda), y = expression(hat(pi)[0](lambda)), title = title) +
    ggplot2::ylim(0, 1)
}

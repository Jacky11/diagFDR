#' Storey pi0(lambda) tail-uniformity diagnostic
#'
#' Estimates pi0(lambda) = #{p > lambda} / ((1-lambda)*m) across a grid of lambdas.
#' Intended as a plausibility diagnostic for (pseudo-)p-values.
#'
#' @param x A dfdr_tbl (or data.frame) containing column \code{p}.
#' @param lambdas Numeric vector in [0,1). Typical range: 0.5--0.95.
#' @param stratify Optional character vector of column names for stratification.
#' @param min_n Minimum number of finite p-values required per stratum.
#' @param clamp Logical; if TRUE clamp pi0 into `[0,1]` for reporting (default TRUE).
#'
#' @return A list with elements \code{pi0} (tibble) and \code{summary} (tibble).
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

#' Plot Storey pi0(lambda) curve
#'
#' @param pi0_tbl Output \code{$pi0} from \code{\link{dfdr_pi0_storey}}.
#' @param title Plot title.
#' @return A ggplot object.
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

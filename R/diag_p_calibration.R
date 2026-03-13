#' P-value calibration diagnostic (ECDF vs uniform; stratified)
#'
#' Computes empirical CDF values Fhat(u) for a grid of u and reports
#' inflation metrics such as max(Fhat(u) - u) in a decision-relevant region.
#' Intended as a plausibility diagnostic for (pseudo-)p-values.
#'
#' @param x A dfdr_tbl (or data.frame) containing columns \code{id} and \code{p}.
#' @param u_grid Numeric vector in (0,1]. Grid of u values for evaluating the ECDF.
#' @param u_max Numeric in (0,1]. Upper bound of u used for inflation summaries (default 0.1).
#' @param stratify Optional character vector of column names used to stratify diagnostics.
#' @param min_n Minimum number of finite p-values required per stratum.
#'
#' @return A list with elements \code{ecdf} (tibble) and \code{summary} (tibble).
#' @export
dfdr_p_calibration <- function(x,
                               u_grid = c(seq(0.001, 0.1, by = 0.001), seq(0.11, 1, by = 0.01)),
                               u_max = 0.1,
                               stratify = NULL,
                               min_n = 200) {
  x <- .safe_df(x)
  .check_has_cols(x, c("id", "p"))

  if (!is.numeric(u_grid) || any(!is.finite(u_grid)) || any(u_grid <= 0) || any(u_grid > 1)) {
    rlang::abort("`u_grid` must be numeric values in (0,1].")
  }
  u_grid <- sort(unique(u_grid))
  if (!is.numeric(u_max) || length(u_max) != 1 || !is.finite(u_max) || u_max <= 0 || u_max > 1) {
    rlang::abort("`u_max` must be a single number in (0,1].")
  }

  strata_df <- .get_strata_df(x, stratify)
  if (is.null(strata_df)) {
    x$stratum <- "all"
  } else {
    x$stratum <- .make_stratum_id(strata_df)
  }

  out_ecdf <- list()
  out_sum <- list()

  for (s in unique(x$stratum)) {
    xs <- x[x$stratum == s, , drop = FALSE]
    p <- xs$p
    p <- p[is.finite(p)]
    n <- length(p)

    if (n < min_n) {
      out_sum[[s]] <- tibble::tibble(
        stratum = s, n = n,
        max_inflation = NA_real_,
        auc_inflation = NA_real_,
        note = paste0("Skipped: n < ", min_n)
      )
      next
    }

    Fhat <- vapply(u_grid, function(u) mean(p <= u), numeric(1))

    ecdf_tbl <- tibble::tibble(
      stratum = s,
      u = u_grid,
      Fhat = Fhat,
      n = n
    )
    out_ecdf[[s]] <- ecdf_tbl

    sel <- u_grid <= u_max
    infl <- (Fhat - u_grid)
    max_infl <- max(infl[sel], na.rm = TRUE)

    # Positive inflation AUC on [0, u_max], normalized
    u_sel <- u_grid[sel]
    infl_pos <- pmax(infl[sel], 0)
    if (length(u_sel) >= 2) {
      auc <- sum(diff(u_sel) * (head(infl_pos, -1) + tail(infl_pos, -1)) / 2) / max(u_sel)
    } else {
      auc <- NA_real_
    }

    out_sum[[s]] <- tibble::tibble(
      stratum = s,
      n = n,
      max_inflation = max_infl,
      auc_inflation = auc,
      note = NA_character_
    )
  }

  ecdf_all <- if (length(out_ecdf)) dplyr::bind_rows(out_ecdf) else tibble::tibble()
  sum_all  <- dplyr::bind_rows(out_sum)

  list(ecdf = ecdf_all, summary = sum_all)
}

#' Plot p-value calibration (ECDF minus uniform)
#'
#' @param ecdf_tbl Output \code{$ecdf} from \code{\link{dfdr_p_calibration}}.
#' @param title Plot title.
#' @return A ggplot object.
#' @export
dfdr_plot_p_calibration <- function(ecdf_tbl,
                                    title = "P-value calibration: ECDF(p) - u") {
  .check_has_cols(ecdf_tbl, c("stratum", "u", "Fhat"))
  ggplot2::ggplot(ecdf_tbl, ggplot2::aes(x = .data$u, y = .data$Fhat - .data$u)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ stratum, scales = "free_y") +
    ggplot2::labs(x = "u", y = "Fhat(u) - u", title = title)
}


#' P-value calibration plot (cp4p-style), with multiple pi0 reference curves
#'
#' Plots ECDF of (1-p) against reference curves derived from cp4p::estim.pi0()
#' using multiple pi0 estimation methods.
#'
#' @param x A dfdr_tbl (or data.frame) containing columns id, is_decoy, p.
#' @param sel Subset: "all", "decoy", or "target".
#' @param nbins,pz Passed to cp4p::estim.pi0().
#' @param step Grid step for abs = 1-p.
#' @param return_data If TRUE, return plot + computed data.
#' @return ggplot (or list if return_data=TRUE)
#' @export
dfdr_plot_p_calibration2 <- function(x,
                                     sel = c("all", "decoy", "target"),
                                     nbins = 20,
                                     pz = 0.05,
                                     step = 0.001,
                                     return_data = FALSE) {
  x <- .safe_df(x)
  .check_has_cols(x, c("id", "is_decoy", "p"))

  sel <- match.arg(tolower(sel), c("all", "decoy", "target"))

  xx <- x
  if (sel == "decoy")  xx <- dplyr::filter(xx, .data$is_decoy)
  if (sel == "target") xx <- dplyr::filter(xx, !.data$is_decoy)

  p <- suppressWarnings(as.numeric(xx$p))
  p <- p[is.finite(p) & p >= 0 & p <= 1]

  if (length(p) < 10) rlang::abort("Too few valid p-values to plot calibration.")
  if (all(p == 0) || all(p == 1)) rlang::abort("Degenerate p-values (all 0 or all 1); cannot plot calibration.")

  cmax <- max(p)
  if (!is.finite(cmax) || cmax <= 0) rlang::abort("max(p) must be > 0 for cp4p-style reference curves.")

  one_minus_p <- 1 - p
  Fr <- stats::ecdf(one_minus_p)

  minabs <- max(0, min(one_minus_p) - 0.05)
  abs_grid <- seq(minabs, 1, by = step)
  fc <- Fr(abs_grid)

  df_main <- tibble::tibble(abs = abs_grid, fc = fc)

  # pi0 estimates (ALL methods)
  r <- cp4p::estim.pi0(p = p, pi0.method = "ALL", nbins = nbins, pz = pz)
  pi0_vec <- as.numeric(r$pi0)

  methods <- c("st.spline","st.boot","jiang","histo","langaas","pounds","abh","slim")
  if (length(pi0_vec) != length(methods)) {
    rlang::abort("cp4p::estim.pi0(pi0.method='ALL') did not return 8 estimates as expected.")
  }

  # dr(abs) = pi0 * (abs - 1 + c) / c, truncated at 0
  df_ref <- purrr::map_dfr(seq_along(methods), function(i) {
    dr_vec <- pi0_vec[i] * (abs_grid - 1 + cmax) / cmax
    dr_vec <- pmax(dr_vec, 0)
    tibble::tibble(
      abs = abs_grid,
      dr = dr_vec,
      method = methods[i],
      pi0 = pi0_vec[i]
    )
  })

  pal <- c(
    "st.spline" = "chartreuse2",
    "st.boot"   = "chartreuse4",
    "jiang"     = "darkgreen",
    "histo"     = "darkorange",
    "langaas"   = "darkorange3",
    "pounds"    = "red3",
    "abh"       = "deepskyblue3",
    "slim"      = "blue3"
  )

  g <- ggplot2::ggplot(df_main, ggplot2::aes(x = .data$abs, y = .data$fc)) +
    ggplot2::geom_line(linewidth = 0.9, colour = "black") +
    ggplot2::geom_line(
      data = df_ref,
      ggplot2::aes(x = .data$abs, y = .data$dr, colour = .data$method, linetype = .data$method),
      linewidth = 0.8,
      alpha = 0.9
    ) +
    ggplot2::scale_colour_manual(values = pal, drop = FALSE) +
    ggplot2::coord_cartesian(ylim = c(0, 1.05)) +
    ggplot2::labs(
      title = "Calibration plot (cp4p-style): ECDF of (1-p) with pi0 reference curves",
      subtitle = paste0("Subset: ", sel, " | n=", length(p)),
      x = "1 - p.value",
      y = "Cumulative Distribution Function of 1 - p.value",
      colour = "pi0 method",
      linetype = "pi0 method"
    ) +
    ggplot2::theme_minimal()

  if (isTRUE(return_data)) {
    return(list(plot = g, pi0 = stats::setNames(pi0_vec, methods), data = list(main = df_main, ref = df_ref)))
  }

  g
}

#' Plot decoy support D_alpha versus alpha
#'
#' @param stab_tbl A stability table (e.g. from [dfdr_curve_stability()]) containing
#'   columns `alpha`, `D_alpha`, and `list`.
#' @param xlab X-axis label.
#' @param title Plot title.
#'
#' @return A ggplot object.
#' @export
dfdr_plot_dalpha <- function(stab_tbl, xlab = "alpha (log10)", title = "Decoy support D_alpha vs alpha") {
  ggplot2::ggplot(stab_tbl, ggplot2::aes(x = alpha, y = D_alpha, color = .data$list)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = title, x = xlab, y = "D_alpha", color = "")
}

#' Plot CV_hat versus alpha
#'
#' @param stab_tbl A stability table (e.g. from \code{\link{dfdr_curve_stability}}) with
#'   columns \code{alpha}, \code{CV_hat}, and \code{list}.
#' @param xlab X-axis label.
#' @param title Plot title.
#'
#' @return A ggplot object.
#' @export
dfdr_plot_cv <- function(stab_tbl, xlab = "alpha (log10)", title = "Stability proxy CV_hat vs alpha") {
  ggplot2::ggplot(stab_tbl, ggplot2::aes(x = alpha, y = CV_hat, color = .data$list)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = title, x = xlab, y = "CV_hat", color = "")
}

#' Plot local decoy support D_alpha_win versus alpha
#'
#' @param dwin_tbl A local-tail table (e.g. from [dfdr_local_tail()]) containing
#'   columns `alpha`, `D_alpha_win`, and `list`.
#' @param xlab X-axis label.
#' @param title Plot title.
#'
#' @return A ggplot object.
#' @export
dfdr_plot_dwin <- function(dwin_tbl, xlab = "alpha (log10)", title = "Local decoy support D_alpha_win vs alpha") {
  ggplot2::ggplot(dwin_tbl, ggplot2::aes(x = alpha, y = D_alpha_win, color = .data$list)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = title, x = xlab, y = "D_alpha_win", color = "")
}

#' Plot threshold elasticity (Jaccard) versus alpha
#'
#' @param el_tbl An elasticity table (e.g. from [dfdr_elasticity()]) containing
#'   columns `alpha`, `jaccard`, and `list`.
#' @param xlab X-axis label.
#' @param title Plot title.
#'
#' @return A ggplot object.
#' @export
dfdr_plot_elasticity <- function(el_tbl, xlab = "alpha (log10)", title = "Threshold elasticity (Jaccard) vs alpha") {
  ggplot2::ggplot(el_tbl, ggplot2::aes(x = alpha, y = jaccard, color = .data$list)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = title, x = xlab, y = "Jaccard", color = "")
}

#' Plot PEP reliability
#'
#' @param bins_tbl A PEP reliability bins table (e.g. the `$bins` output of
#'   [dfdr_pep_reliability()]) containing columns `pep_mean`, `decoy_rate`, and `n`.
#' @param title Plot title.
#'
#' @return A ggplot object.
#' @export
dfdr_plot_pep_reliability <- function(bins_tbl, title = "PEP reliability") {
  ggplot2::ggplot(bins_tbl, ggplot2::aes(x = pep_mean, y = decoy_rate, size = n)) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = title, x = "Mean predicted PEP", y = "Observed decoy fraction", size = "n")
}

#' Plot equal-chance plausibility by q-value bands
#'
#' @param bands_tbl A q-band table (e.g. the `$bands` output of
#'   [dfdr_equal_chance_qbands()]) containing columns `q_mean`, `decoy_frac`, and `n`.
#' @param title Plot title.
#'
#' @return A ggplot object.
#' @export
dfdr_plot_equal_chance <- function(bands_tbl, title = "Equal-chance plausibility: decoy fraction by q-band") {
  ggplot2::ggplot(bands_tbl, ggplot2::aes(x = q_mean, y = decoy_frac, size = n)) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::geom_line(alpha = 0.5) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = title, x = "Mean q-value in band", y = "Decoy fraction", size = "n")
}


#' Plot target vs decoy score distributions
#'
#' @param x A dfdr_tbl object
#' @param title Plot title
#' @return A ggplot object with target distribution in blue and decoy distribution in red
#' @export
dfdr_plot_score_distributions <- function(x, title = "Target vs Decoy Score Distributions") {
  validate_dfdr_tbl(x)

  # Determine which variable to plot
  use_score <- "score" %in% names(x) && !all(is.na(x$score))
  use_q <- !use_score  # Fallback to q-values if no score

  if (use_score) {
    plot_data <- x |>
      dplyr::filter(!is.na(score)) |>
      dplyr::mutate(
        type = ifelse(is_decoy, "Decoy", "Target"),
        value = score
      )
    x_label <- "Score"
    plot_title <- title
    use_log <- FALSE
  } else {
    plot_data <- x |>
      dplyr::filter(!is.na(q), q > 0) |>  # Filter out q=0 for log scale
      dplyr::mutate(
        type = ifelse(is_decoy, "Decoy", "Target"),
        value = q
      )
    x_label <- "q-value"
    plot_title <- sub("Score [Dd]istributions?", "q-value distributions", title)
    use_log <- TRUE
  }

  if (nrow(plot_data) == 0) {
    rlang::warn("No valid data for distribution plot; returning empty plot.")
    return(ggplot2::ggplot() + ggplot2::labs(title = paste(plot_title, "(no data)")))
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value, fill = type, color = type)) +
    ggplot2::geom_density(alpha = 0.3, linewidth = 0.8) +
    ggplot2::scale_fill_manual(values = c("Target" = "#2E86AB", "Decoy" = "#A23B72")) +
    ggplot2::scale_color_manual(values = c("Target" = "#2E86AB", "Decoy" = "#A23B72")) +
    ggplot2::labs(
      title = plot_title,
      x = x_label,
      y = "Density",
      fill = "",
      color = ""
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "top",
      plot.title = ggplot2::element_text(face = "bold")
    )

  # Use log scale for q-values
  if (use_log) {
    p <- p +
      ggplot2::scale_x_log10(
        labels = scales::label_number(accuracy = 0.001)
      ) +
      ggplot2::annotation_logticks(sides = "b")
  }

  p
}


#' Plot scope disagreement as Jaccard overlap heatmap
#'
#' @param xs Named list of dfdr_tbl objects
#' @param alpha Threshold for accepted targets
#' @return A ggplot heatmap
#' @export
dfdr_plot_scope_disagreement_matrix <- function(xs, alpha = 0.01) {

  # Helper to extract base precursor ID (strip run prefix if present)
  extract_base_id <- function(ids) {
    # If IDs contain "||", extract part after it (runxprecursor format)
    # Otherwise use as-is (precursor format)
    ifelse(grepl("\\|\\|", ids),
           sub("^.*\\|\\|", "", ids),  # Everything after ||
           ids)                         # Keep as-is
  }

  # Get accepted target sets with normalized IDs
  accepted <- purrr::map(xs, ~{
    .x |>
      dplyr::filter(!is_decoy, q <= alpha) |>
      dplyr::pull(id) |>
      extract_base_id() |>
      unique()
  })

  # Compute pairwise Jaccard
  n_lists <- length(accepted)
  list_names <- names(accepted)

  jaccard_mat <- matrix(NA, n_lists, n_lists, dimnames = list(list_names, list_names))

  for (i in seq_len(n_lists)) {
    for (j in seq_len(n_lists)) {
      if (i == j) {
        jaccard_mat[i, j] <- 1.0
      } else {
        intersection <- length(intersect(accepted[[i]], accepted[[j]]))
        union_size <- length(union(accepted[[i]], accepted[[j]]))
        jaccard_mat[i, j] <- if (union_size > 0) intersection / union_size else 0
      }
    }
  }

  # Convert to long format for ggplot
  jaccard_df <- as.data.frame(jaccard_mat) |>
    tibble::rownames_to_column("list1") |>
    tidyr::pivot_longer(-list1, names_to = "list2", values_to = "jaccard")

  ggplot2::ggplot(jaccard_df, ggplot2::aes(x = list1, y = list2, fill = jaccard)) +
    ggplot2::geom_tile(color = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", jaccard)),
                       size = 5, fontface = "bold") +
    ggplot2::scale_fill_gradient2(
      low = "#d73027", mid = "#fee08b", high = "#1a9850",
      midpoint = 0.85, limits = c(0, 1),
      name = "Jaccard\noverlap"
    ) +
    ggplot2::labs(
      title = sprintf("Scope disagreement at α = %.2g", alpha),
      subtitle = "Jaccard overlap of accepted target sets (precursor-level comparison)",
      x = "", y = ""
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = ggplot2::element_text(size = 11),
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 11)
    ) +
    ggplot2::coord_fixed()
}


#' Plot PEP density by decoy/target
#'
#' @param x A dfdr_tbl with pep.
#' @param pep_max Max PEP displayed (default 0.5).
#' @param bw Bandwidth passed to stats::density() via geom_density (optional).
#' @param adjust Smoothing adjustment passed to geom_density (default 1).
#' @param title Plot title.
#' @return ggplot
#' @export
dfdr_plot_pep_density_by_decoy <- function(x, pep_max = 0.5,
                                           bw = "nrd0",
                                           adjust = 1,
                                           title = "PEP distribution (density): targets vs decoys") {
  validate_dfdr_tbl(x)
  if (!"pep" %in% names(x)) rlang::abort("`x` must contain a `pep` column.")

  df <- x |>
    dplyr::filter(is.finite(pep), pep >= 0, pep <= pep_max) |>
    dplyr::mutate(class = ifelse(is_decoy, "Decoy", "Target"))

  ggplot2::ggplot(df, ggplot2::aes(x = pep, color = class, fill = class)) +
    ggplot2::geom_density(
      linewidth = 1,
      alpha = 0.25,
      bw = bw,
      adjust = adjust,
      na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(values = c(Target = "#2b8cbe", Decoy = "#d7301f")) +
    ggplot2::scale_fill_manual(values = c(Target = "#2b8cbe", Decoy = "#d7301f")) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("Density truncated to pep ≤ %.2f", pep_max),
      x = "PEP", y = "Density",
      color = NULL, fill = NULL
    ) +
    ggplot2::theme_minimal()
}

#' Plot target-focused PEP reliability (TDC-style)
#'
#' @param rel_tbl Output of dfdr_pep_reliability_tdc()
#' @return ggplot
#' @export
dfdr_plot_pep_reliability_tdc <- function(rel_tbl,
                                          title = "PEP reliability (targets; TDC-style D/T error proxy)") {
  ggplot2::ggplot(rel_tbl, ggplot2::aes(x = pep_mean_targets, y = err_hat_target)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(ggplot2::aes(size = n_target), alpha = 0.8) +
    ggplot2::scale_size_continuous(range = c(2, 10)) +
    ggplot2::coord_cartesian(xlim = c(0, max(rel_tbl$pep_mean_targets, na.rm = TRUE)),
                             ylim = c(0, 1)) +
    ggplot2::labs(
      title = title,
      x = "Mean predicted PEP (targets in bin)",
      y = "Estimated false fraction among targets (min(1, (D+add)/T))",
      size = "Targets/bin"
    ) +
    ggplot2::theme_minimal()
}


#' Plot (pseudo-)p-value density by decoy/target
#'
#' @param x A dfdr_tbl with column p and is_decoy.
#' @param p_max Max p-value displayed (default 1).
#' @param bw Bandwidth passed to geom_density (optional).
#' @param adjust Smoothing adjustment (default 1).
#' @param title Plot title.
#' @param log10_x If TRUE, plot -log10(p) instead of p (often more informative).
#' @param eps Lower bound to avoid Inf when log10_x=TRUE.
#' @return ggplot
#' @export
dfdr_plot_p_density_by_decoy <- function(x,
                                         p_max = 1,
                                         bw = "nrd0",
                                         adjust = 1,
                                         title = "(Pseudo-)p-value distribution (density): targets vs decoys",
                                         log10_x = FALSE,
                                         eps = 1e-300) {
  validate_dfdr_tbl(x)
  if (!"p" %in% names(x)) rlang::abort("`x` must contain a `p` column.")

  df <- x |>
    dplyr::filter(is.finite(p), p > 0, p <= p_max) |>
    dplyr::mutate(
      class = ifelse(is_decoy, "Decoy", "Target"),
      xval = if (log10_x) -log10(pmax(p, eps)) else p
    )

  xlab <- if (log10_x) expression(-log[10](p)) else "p"

  ggplot2::ggplot(df, ggplot2::aes(x = xval, color = class, fill = class)) +
    ggplot2::geom_density(
      linewidth = 1,
      alpha = 0.25,
      bw = bw,
      adjust = adjust,
      na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(values = c(Target = "#2b8cbe", Decoy = "#d7301f")) +
    ggplot2::scale_fill_manual(values = c(Target = "#2b8cbe", Decoy = "#d7301f")) +
    ggplot2::labs(
      title = title,
      subtitle = if (log10_x) sprintf("Using -log10(p); computed p truncated to p ≤ %.2g", p_max)
      else sprintf("Computed p truncated to p ≤ %.2g", p_max),
      x = xlab, y = "Density",
      color = NULL, fill = NULL
    ) +
    ggplot2::theme_minimal()
}




#' Plot Storey-style estimated FDR curve: pi0 * m * t / R(t)
#'
#' Given p-values (or pseudo-p-values), estimates pi0 using cp4p::estim.pi0()
#' and plots the curve:
#'   FDR_hat(t) = (pi0_hat * m * t) / R(t),
#' where R(t) = #{p <= t}, m = total number of tested hypotheses (after filtering).
#'
#' @param x A dfdr_tbl (or data.frame) containing columns \code{id}, \code{p}, \code{is_decoy}.
#' @param sel Subset to use: "all", "target", or "decoy" (case-insensitive).
#'   For "target", keeps !is_decoy; for "decoy", keeps is_decoy.
#' @param pi0.method Passed to cp4p::estim.pi0(). One of "st.spline","st.boot","jiang",
#'   "histo","langaas","pounds","abh","slim". (Do not use "ALL" here.)
#' @param nbins,pz Passed to cp4p::estim.pi0().
#' @param t_grid Grid of t values in (0,1] at which to evaluate the curve.
#' @param add_one If TRUE, uses R(t) := #{p <= t} + 1 to avoid division by zero at tiny t.
#' @param cap_one If TRUE, caps FDR_hat at 1.
#'
#' @return A list with:
#'   \itemize{
#'     \item plot: ggplot object
#'     \item data: tibble with t, R, fdr_hat
#'     \item pi0_hat: numeric
#'     \item m: integer
#'   }
#' @export
dfdr_plot_fdrhat_pi0 <- function(x,
                                 sel = c("all","target","decoy"),
                                 pi0.method = "pounds",
                                 nbins = 20,
                                 pz = 0.05,
                                 t_grid = 10^seq(-6, 0, length.out = 400),
                                 r_min = 1,
                                 cap_one = TRUE) {
  x <- .safe_df(x)
  .check_has_cols(x, c("id","p","is_decoy"))

  sel <- match.arg(tolower(sel), c("all","target","decoy"))
  xx <- x
  if (sel == "target") xx <- dplyr::filter(xx, !.data$is_decoy)
  if (sel == "decoy")  xx <- dplyr::filter(xx,  .data$is_decoy)

  p <- suppressWarnings(as.numeric(xx$p))
  p <- p[is.finite(p) & p >= 0 & p <= 1]
  m <- length(p)
  if (m < 50) rlang::abort("Too few valid p-values to estimate pi0/FDR curve reliably.")

  r <- cp4p::estim.pi0(p = p, pi0.method = pi0.method, nbins = nbins, pz = pz)
  pi0_hat <- as.numeric(r$pi0)

  p_sorted <- sort(p)
  t_grid <- t_grid[is.finite(t_grid) & t_grid > 0 & t_grid <= 1]
  R <- findInterval(t_grid, p_sorted, rightmost.closed = TRUE)

  fdr_hat <- rep(NA_real_, length(t_grid))
  ok <- R >= r_min
  fdr_hat[ok] <- (pi0_hat * m * t_grid[ok]) / R[ok]
  if (cap_one) fdr_hat <- pmin(1, fdr_hat)

  df <- tibble::tibble(t = t_grid, R = R, fdr_hat = fdr_hat)

  g <- ggplot2::ggplot(df, ggplot2::aes(t, fdr_hat)) +
    ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
    ggplot2::scale_x_log10() +
    ggplot2::coord_cartesian(ylim = c(0,1)) +
    ggplot2::labs(
      title = expression(hat(FDR)(t) == frac(hat(pi)[0] * m * t, R(t))),
      subtitle = paste0("subset=", sel, " | pi0_hat=", signif(pi0_hat,3),
                        " | m=", m, " | r_min=", r_min),
      x = "t (p-value threshold, log scale)",
      y = expression(hat(FDR)(t))
    ) +
    ggplot2::theme_minimal()

  list(plot = g, data = df, pi0_hat = pi0_hat, m = m)
}


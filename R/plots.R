#' Plot decoy support \eqn{D_\alpha} versus \code{alpha}
#'
#' Plots the number of accepted decoys at each threshold \code{alpha}. This is a
#' key stability indicator: very small \eqn{D_\alpha} implies granular/unstable
#' behaviour of target-decoy based estimates.
#'
#' @param stab_tbl A stability table (e.g. from \code{\link{dfdr_curve_stability}})
#'   containing columns \code{alpha}, \code{D_alpha}, and \code{list}.
#' @param xlab Character. X-axis label.
#' @param title Character. Plot title.
#'
#' @return A \link[ggplot2:ggplot]{ggplot} object.
#'
#' @examples
#' library(tibble)
#'
#' stab_tbl <- tibble(
#'   list = rep(c("A", "B"), each = 4),
#'   alpha = rep(c(1e-3, 2e-3, 5e-3, 1e-2), times = 2),
#'   D_alpha = c(5, 8, 20, 40, 2, 3, 6, 10)
#' )
#' dfdr_plot_dalpha(stab_tbl)
#'
#' @export
dfdr_plot_dalpha <- function(stab_tbl, xlab = "alpha (log10)", title = "Decoy support D_alpha vs alpha") {
  ggplot2::ggplot(stab_tbl, ggplot2::aes(x = alpha, y = D_alpha, color = .data$list)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = title, x = xlab, y = "D_alpha", color = "")
}

#' Plot \code{CV_hat} versus \code{alpha}
#'
#' Plots a stability proxy (\code{CV_hat}) against the threshold \code{alpha} for
#' one or more lists.
#'
#' @param stab_tbl A stability table (e.g. from \code{\link{dfdr_curve_stability}})
#'   with columns \code{alpha}, \code{CV_hat}, and \code{list}.
#' @param xlab Character. X-axis label.
#' @param title Character. Plot title.
#'
#' @return A \link[ggplot2:ggplot]{ggplot} object.
#'
#' @examples
#' library(tibble)
#'
#' stab_tbl <- tibble(
#'   list = rep(c("A", "B"), each = 4),
#'   alpha = rep(c(1e-3, 2e-3, 5e-3, 1e-2), times = 2),
#'   CV_hat = c(0.30, 0.22, 0.15, 0.12, 0.40, 0.28, 0.20, 0.18)
#' )
#' dfdr_plot_cv(stab_tbl)
#'
#' @export
dfdr_plot_cv <- function(stab_tbl, xlab = "alpha (log10)", title = "Stability proxy CV_hat vs alpha") {
  ggplot2::ggplot(stab_tbl, ggplot2::aes(x = alpha, y = CV_hat, color = .data$list)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = title, x = xlab, y = "CV_hat", color = "")
}

#' Plot local decoy support \code{D_alpha_win} versus \code{alpha}
#'
#' Plots the number of decoys in a right-neighborhood above each threshold
#' (as computed by \code{\link{dfdr_local_tail}}). This approximates how well
#' supported the boundary is by nearby decoys.
#'
#' @param dwin_tbl A local-tail table (e.g. from \code{\link{dfdr_local_tail}})
#'   containing columns \code{alpha}, \code{D_alpha_win}, and \code{list}.
#' @param xlab Character. X-axis label.
#' @param title Character. Plot title.
#'
#' @return A \link[ggplot2:ggplot]{ggplot} object.
#'
#' @examples
#' library(tibble)
#'
#' dwin_tbl <- tibble(
#'   list = rep(c("A", "B"), each = 4),
#'   alpha = rep(c(1e-3, 2e-3, 5e-3, 1e-2), times = 2),
#'   D_alpha_win = c(1, 2, 6, 15, 0, 1, 2, 4)
#' )
#' dfdr_plot_dwin(dwin_tbl)
#'
#' @export
dfdr_plot_dwin <- function(dwin_tbl, xlab = "alpha (log10)", title = "Local decoy support D_alpha_win vs alpha") {
  ggplot2::ggplot(dwin_tbl, ggplot2::aes(x = alpha, y = D_alpha_win, color = .data$list)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = title, x = xlab, y = "D_alpha_win", color = "")
}

#' Plot threshold elasticity (Jaccard) versus \code{alpha}
#'
#' Plots Jaccard overlap values (as returned by \code{\link{dfdr_elasticity}})
#' against \code{alpha}. Lower overlap indicates greater sensitivity of the
#' accepted set to small changes of the threshold.
#'
#' @param el_tbl An elasticity table (e.g. from \code{\link{dfdr_elasticity}})
#'   containing columns \code{alpha}, \code{jaccard}, and \code{list}.
#' @param xlab Character. X-axis label.
#' @param title Character. Plot title.
#'
#' @return A \link[ggplot2:ggplot]{ggplot} object.
#'
#' @examples
#' library(tibble)
#'
#' el_tbl <- tibble(
#'   list = rep(c("A", "B"), each = 4),
#'   alpha = rep(c(1e-3, 2e-3, 5e-3, 1e-2), times = 2),
#'   jaccard = c(0.95, 0.93, 0.90, 0.88, 0.92, 0.90, 0.86, 0.82)
#' )
#' dfdr_plot_elasticity(el_tbl)
#'
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
#' Visualises PEP reliability by plotting the observed decoy fraction against the
#' mean predicted PEP in each bin (from \code{\link{dfdr_pep_reliability}}). Point
#' size reflects bin counts.
#'
#' @param bins_tbl A PEP reliability bins table (e.g. \code{out$bins} from
#'   \code{\link{dfdr_pep_reliability}}) containing columns \code{pep_mean},
#'   \code{decoy_rate}, and \code{n}.
#' @param title Character. Plot title.
#'
#' @return A \link[ggplot2:ggplot]{ggplot} object.
#'
#' @examples
#' library(tibble)
#'
#' bins_tbl <- tibble(
#'   pep_mean = c(0.05, 0.15, 0.25, 0.35),
#'   decoy_rate = c(0.06, 0.14, 0.30, 0.40),
#'   n = c(500, 400, 250, 120)
#' )
#' dfdr_plot_pep_reliability(bins_tbl)
#'
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
#' Visualises decoy fractions by q-value band (as produced by
#' \code{\link{dfdr_equal_chance_qbands}}). A horizontal reference at 0.5 indicates
#' the expected decoy fraction under an "equal-chance" region assumption.
#'
#' @param bands_tbl A q-band table (e.g. \code{out$bands} from
#'   \code{\link{dfdr_equal_chance_qbands}}) containing columns \code{q_mean},
#'   \code{decoy_frac}, and \code{n}.
#' @param title Character. Plot title.
#'
#' @return A \link[ggplot2:ggplot]{ggplot} object.
#'
#' @examples
#' library(tibble)
#'
#' bands_tbl <- tibble(
#'   q_mean = c(0.05, 0.15, 0.30, 0.45),
#'   decoy_frac = c(0.10, 0.20, 0.45, 0.52),
#'   n = c(2000, 1500, 800, 400)
#' )
#' dfdr_plot_equal_chance(bands_tbl)
#'
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
#' Plots density distributions of scores (if available) or q-values (fallback) for
#' targets vs decoys. This provides a quick sanity check that decoys are shifted
#' toward lower scores (or higher q-values).
#'
#' @param x A \code{dfdr_tbl}.
#' @param title Character. Plot title.
#'
#' @return A \link[ggplot2:ggplot]{ggplot} object.
#'
#' @examples
#' library(tibble)
#'
#' set.seed(1)
#' n <- 4000
#' df <- tibble(
#'   id = as.character(seq_len(n)),
#'   run = "run1",
#'   is_decoy = sample(c(FALSE, TRUE), n, replace = TRUE, prob = c(0.95, 0.05)),
#'   score = rnorm(n),
#'   q = pmin(1, rank(-score) / n),
#'   pep = NA_real_
#' )
#' x <- as_dfdr_tbl(df, unit = "psm", scope = "global", q_source = "toy")
#' dfdr_plot_score_distributions(x)
#'
#' @export
dfdr_plot_score_distributions <- function(x, title = "Target vs Decoy Score Distributions") {
  validate_dfdr_tbl(x)

  use_score <- "score" %in% names(x) && !all(is.na(x$score))

  if (use_score) {
    plot_data <- x |>
      dplyr::filter(is.finite(.data$score)) |>
      dplyr::mutate(
        type = dplyr::if_else(.data$is_decoy, "Decoy", "Target"),
        value = .data$score
      )
    x_label <- "Score"
    plot_title <- title
    use_log <- FALSE
  } else {
    plot_data <- x |>
      dplyr::filter(is.finite(.data$q), .data$q > 0) |>
      dplyr::mutate(
        type = dplyr::if_else(.data$is_decoy, "Decoy", "Target"),
        value = .data$q
      )
    x_label <- "q-value"
    plot_title <- sub("Score [Dd]istributions?", "q-value distributions", title)
    use_log <- TRUE
  }

  if (nrow(plot_data) == 0) {
    rlang::warn("No valid data for distribution plot; returning empty plot.")
    return(ggplot2::ggplot() + ggplot2::labs(title = paste(plot_title, "(no data)")))
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$value, fill = .data$type, color = .data$type)) +
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

  if (use_log) {
    p <- p +
      ggplot2::scale_x_log10(labels = scales::label_number(accuracy = 0.001)) +
      ggplot2::annotation_logticks(sides = "b")
  }

  p
}

#' Plot scope disagreement as a Jaccard overlap heatmap
#'
#' Computes accepted target ID sets for each list at a fixed \code{alpha} and
#' visualises pairwise Jaccard overlaps as a heatmap. IDs are compared at the
#' "base" level by stripping any \code{"run||"} prefix.
#'
#' @param xs Named list of \code{dfdr_tbl} objects.
#' @param alpha Numeric threshold used to define accepted targets (\code{q <= alpha}).
#'
#' @return A \link[ggplot2:ggplot]{ggplot} object (heatmap).
#'
#' @examples
#' library(tibble)
#'
#' set.seed(1)
#' n <- 2000
#'
#' make_x <- function(label) {
#'   df <- tibble(
#'     id = paste0("run1||", as.character(seq_len(n))),
#'     run = "run1",
#'     is_decoy = sample(c(FALSE, TRUE), n, replace = TRUE, prob = c(0.95, 0.05)),
#'     score = rnorm(n),
#'     q = stats::runif(n, 0, 0.05),
#'     pep = NA_real_
#'   )
#'   as_dfdr_tbl(df, unit = "psm", scope = label, q_source = "toy")
#' }
#'
#' xs <- list(A = make_x("A"), B = make_x("B"), C = make_x("C"))
#' dfdr_plot_scope_disagreement_matrix(xs, alpha = 0.01)
#'
#' @export
dfdr_plot_scope_disagreement_matrix <- function(xs, alpha = 0.01) {
  extract_base_id <- function(ids) {
    ifelse(grepl("\\|\\|", ids), sub("^.*\\|\\|", "", ids), ids)
  }

  accepted <- purrr::map(xs, ~{
    .x |>
      dplyr::filter(!.data$is_decoy, .data$q <= alpha) |>
      dplyr::pull(.data$id) |>
      extract_base_id() |>
      unique()
  })

  n_lists <- length(accepted)
  list_names <- names(accepted)

  jaccard_mat <- matrix(NA_real_, n_lists, n_lists, dimnames = list(list_names, list_names))

  for (i in seq_len(n_lists)) {
    for (j in seq_len(n_lists)) {
      if (i == j) {
        jaccard_mat[i, j] <- 1.0
      } else {
        inter <- length(intersect(accepted[[i]], accepted[[j]]))
        uni <- length(union(accepted[[i]], accepted[[j]]))
        jaccard_mat[i, j] <- if (uni > 0) inter / uni else 0
      }
    }
  }

  jaccard_df <- as.data.frame(jaccard_mat) |>
    tibble::rownames_to_column("list1") |>
    tidyr::pivot_longer(
      cols = -dplyr::all_of("list1"),
      names_to = "list2",
      values_to = "jaccard"
    )

  ggplot2::ggplot(jaccard_df, ggplot2::aes(x = .data$list1, y = .data$list2, fill = .data$jaccard)) +
    ggplot2::geom_tile(color = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$jaccard)),
                       size = 5, fontface = "bold") +
    ggplot2::scale_fill_gradient2(
      low = "#d73027", mid = "#fee08b", high = "#1a9850",
      midpoint = 0.85, limits = c(0, 1),
      name = "Jaccard\noverlap"
    ) +
    ggplot2::labs(
      title = sprintf("Scope disagreement at alpha = %.2g", alpha),
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
#' Plots density curves of PEP values for targets vs decoys, truncated to
#' \code{pep <= pep_max}.
#'
#' @param x A \code{dfdr_tbl} with a \code{pep} column.
#' @param pep_max Numeric. Maximum PEP displayed (default 0.5).
#' @param bw Bandwidth passed to \code{ggplot2::geom_density()}.
#' @param adjust Numeric smoothing adjustment passed to \code{geom_density()}.
#' @param title Character. Plot title.
#'
#' @return A \link[ggplot2:ggplot]{ggplot} object.
#'
#' @examples
#' library(tibble)
#'
#' set.seed(1)
#' n <- 4000
#' df <- tibble(
#'   id = as.character(seq_len(n)),
#'   run = "run1",
#'   is_decoy = sample(c(FALSE, TRUE), n, replace = TRUE, prob = c(0.95, 0.05)),
#'   score = rnorm(n),
#'   q = stats::runif(n, 0, 0.05),
#'   pep = stats::runif(n)
#' )
#' x <- as_dfdr_tbl(df, unit = "psm", scope = "global", q_source = "toy")
#' dfdr_plot_pep_density_by_decoy(x, pep_max = 0.5)
#'
#' @export
dfdr_plot_pep_density_by_decoy <- function(x, pep_max = 0.5,
                                           bw = "nrd0",
                                           adjust = 1,
                                           title = "PEP distribution (density): targets vs decoys") {
  validate_dfdr_tbl(x)
  if (!"pep" %in% names(x)) rlang::abort("`x` must contain a `pep` column.")

  df <- x |>
    dplyr::filter(is.finite(.data$pep), .data$pep >= 0, .data$pep <= pep_max) |>
    dplyr::mutate(class = dplyr::if_else(.data$is_decoy, "Decoy", "Target"))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$pep, color = .data$class, fill = .data$class)) +
    ggplot2::geom_density(linewidth = 1, alpha = 0.25, bw = bw, adjust = adjust, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(Target = "#2b8cbe", Decoy = "#d7301f")) +
    ggplot2::scale_fill_manual(values = c(Target = "#2b8cbe", Decoy = "#d7301f")) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("Density truncated to pep <= %.2f", pep_max),
      x = "PEP", y = "Density",
      color = NULL, fill = NULL
    ) +
    ggplot2::theme_minimal()
}

#' Plot target-focused PEP reliability (TDC-style)
#'
#' Visualises the output of \code{\link{dfdr_pep_reliability_tdc}} by plotting the
#' estimated target error proxy \code{err_hat_target} against the mean predicted
#' target PEP per bin. Point size reflects the number of targets per bin.
#'
#' @param rel_tbl Output of \code{\link{dfdr_pep_reliability_tdc}}.
#' @param title Character. Plot title.
#'
#' @return A \link[ggplot2:ggplot]{ggplot} object.
#'
#' @examples
#' library(tibble)
#'
#' rel_tbl <- tibble(
#'   bin = c("(0,0.1]", "(0.1,0.2]", "(0.2,0.3]"),
#'   n_target = c(800, 500, 200),
#'   n_decoy = c(20, 30, 40),
#'   pep_mean_targets = c(0.05, 0.15, 0.25),
#'   pep_mean_all = c(0.06, 0.17, 0.28),
#'   err_hat_target = c(0.03, 0.06, 0.20),
#'   pep_bin_mid = c(0.05, 0.15, 0.25)
#' )
#' dfdr_plot_pep_reliability_tdc(rel_tbl)
#'
#' @export
dfdr_plot_pep_reliability_tdc <- function(rel_tbl,
                                          title = "PEP reliability (targets; TDC-style D/T error proxy)") {
  ggplot2::ggplot(rel_tbl, ggplot2::aes(x = .data$pep_mean_targets, y = .data$err_hat_target)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(ggplot2::aes(size = .data$n_target), alpha = 0.8) +
    ggplot2::scale_size_continuous(range = c(2, 10)) +
    ggplot2::coord_cartesian(
      xlim = c(0, max(rel_tbl$pep_mean_targets, na.rm = TRUE)),
      ylim = c(0, 1)
    ) +
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
#' Plots density curves of (pseudo-)p-values for targets vs decoys. Optionally
#' plots \eqn{-\log_{10}(p)} to emphasise small p-values.
#'
#' @param x A \code{dfdr_tbl} with columns \code{p} and \code{is_decoy}.
#' @param p_max Numeric. Maximum p-value displayed (default 1).
#' @param bw Bandwidth passed to \code{ggplot2::geom_density()}.
#' @param adjust Numeric smoothing adjustment passed to \code{geom_density()}.
#' @param title Character. Plot title.
#' @param log10_x Logical. If \code{TRUE}, plot \eqn{-\log_{10}(p)} instead of \code{p}.
#' @param eps Numeric. Lower bound used to avoid \code{Inf} when \code{log10_x = TRUE}.
#'
#' @return A \link[ggplot2:ggplot]{ggplot} object.
#'
#' @examples
#' library(tibble)
#'
#' set.seed(1)
#' n <- 4000
#' df <- tibble(
#'   id = as.character(seq_len(n)),
#'   run = "run1",
#'   is_decoy = sample(c(FALSE, TRUE), n, replace = TRUE, prob = c(0.95, 0.05)),
#'   score = rnorm(n),
#'   q = stats::runif(n, 0, 0.05),
#'   pep = NA_real_,
#'   p = stats::runif(n)
#' )
#' x <- as_dfdr_tbl(df, unit = "psm", scope = "global", q_source = "toy", p_source = "toy")
#' dfdr_plot_p_density_by_decoy(x, log10_x = TRUE)
#'
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
    dplyr::filter(is.finite(.data$p), .data$p > 0, .data$p <= p_max) |>
    dplyr::mutate(
      class = dplyr::if_else(.data$is_decoy, "Decoy", "Target"),
      xval = if (log10_x) -log10(pmax(.data$p, eps)) else .data$p
    )

  xlab <- if (log10_x) expression(-log[10](p)) else "p"

  ggplot2::ggplot(df, ggplot2::aes(x = .data$xval, color = .data$class, fill = .data$class)) +
    ggplot2::geom_density(linewidth = 1, alpha = 0.25, bw = bw, adjust = adjust, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(Target = "#2b8cbe", Decoy = "#d7301f")) +
    ggplot2::scale_fill_manual(values = c(Target = "#2b8cbe", Decoy = "#d7301f")) +
    ggplot2::labs(
      title = title,
      subtitle = if (log10_x) sprintf("Using -log10(p); p truncated to p <= %.2g", p_max)
      else sprintf("p truncated to p <= %.2g", p_max),
      x = xlab, y = "Density",
      color = NULL, fill = NULL
    ) +
    ggplot2::theme_minimal()
}


#' Plot Storey-style estimated FDR curve: \eqn{\hat\pi_0 m t / R(t)}
#'
#' Given p-values (or pseudo-p-values), estimates \eqn{\pi_0} using
#' \code{cp4p::estim.pi0()} and plots
#' \deqn{\widehat{\mathrm{FDR}}(t) = \hat\pi_0 \, m \, t / R(t),}
#' where \eqn{R(t) = \#\{p \le t\}} and \eqn{m} is the number of finite p-values.
#'
#' @param x A \code{dfdr_tbl} (or data.frame) containing columns \code{id},
#'   \code{p}, and \code{is_decoy}.
#' @param sel Character. Subset to use: \code{"all"}, \code{"target"}, or \code{"decoy"}.
#' @param pi0.method Character. Method passed to \code{cp4p::estim.pi0()} (not \code{"ALL"}).
#' @param nbins Integer. Passed to \code{cp4p::estim.pi0()}.
#' @param pz Numeric. Passed to \code{cp4p::estim.pi0()}.
#' @param t_grid Numeric vector of thresholds in \eqn{(0,1]} at which to evaluate the curve.
#' @param r_min Numeric. Minimum value used to determine when to compute the curve
#'   to avoid division by zero at tiny \code{t}.
#' @param cap_one Logical. If \code{TRUE}, caps \code{fdr_hat} at 1.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{plot}{A \link[ggplot2:ggplot]{ggplot} object.}
#'   \item{data}{A \link[tibble:tibble]{tibble} with columns \code{t}, \code{R}, and \code{fdr_hat}.}
#'   \item{pi0_hat}{Numeric scalar \eqn{\hat\pi_0}.}
#'   \item{m}{Integer. Number of finite p-values used.}
#' }
#'
#' @examples
#' library(tibble)
#'
#' if (requireNamespace("cp4p", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 4000
#'   df <- tibble(
#'     id = as.character(seq_len(n)),
#'     is_decoy = sample(c(FALSE, TRUE), n, replace = TRUE, prob = c(0.95, 0.05)),
#'     p = c(stats::runif(3600), stats::rbeta(400, 0.3, 1))
#'   )
#'   out <- dfdr_plot_fdrhat_pi0(df, sel = "all", pi0.method = "pounds", nbins = 10)
#'   out$pi0_hat
#'   out$plot
#' }
#'
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


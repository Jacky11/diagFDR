#' Run a standard set of FDR QC diagnostics
#'
#' Runs stability, calibration, and (optionally) PEP- and p-value-based diagnostics
#' on a named list of \code{dfdr_tbl} objects and returns tables and plots.
#'
#' @param xs Named list of \code{dfdr_tbl} objects.
#' @param alpha_main Headline threshold (e.g. 0.01).
#' @param alphas Numeric vector of thresholds used to build curves.
#' @param eps Relative perturbation used in elasticity (e.g. 0.2).
#' @param win_rel Relative window width for local tail support (e.g. 0.2).
#' @param truncation Truncation policy for local tail windows; see \code{\link{dfdr_local_tail}}.
#' @param k_fixed Fixed +/-K used in sensitivity intervals.
#' @param k_sqrt_mult Multiplier for adaptive +/-ceil(mult*sqrt(D)) intervals.
#' @param qband_breaks Cutpoints for q-band equal-chance diagnostic.
#' @param low_conf Pooled low-confidence region for equal-chance test.
#' @param min_N_equalchance Minimum pooled N for equal-chance test.
#' @param pep_binwidth Bin width for PEP reliability bins.
#' @param pep_n_min Minimum bin size to include PEP bins in IPE summary/plots.
#' @param pep_max Max mean PEP included in IPE summary.
#'
#' @param pep_sanity_thresholds Thresholds for decoy PEP sanity summaries.
#' @param pep_density_max Max PEP displayed in target/decoy density plot.
#' @param pep_rel_tdc_breaks Breaks used for the TDC-style PEP reliability bins.
#' @param pep_rel_tdc_add_decoy Additive correction used in D/T for TDC-style reliability.
#'
#' @param compute_pseudo_pvalues If TRUE, attempt to compute pseudo-p-values from score where p is missing.
#' @param pseudo_pvalue_method Method passed to \code{\link{score_to_pvalue}} (default "decoy_ecdf").
#' @param p_u_grid Grid for p-value ECDF calibration diagnostics.
#' @param p_u_max Upper bound of u used for p-value inflation summaries.
#' @param p_stratify Optional character vector of columns for stratified p-value diagnostics.
#' @param p_min_n_calib Minimum n for p-value calibration per stratum.
#' @param p_lambdas Lambdas for Storey pi0(lambda) diagnostic.
#' @param p_min_n_pi0 Minimum n for pi0 diagnostic per stratum.
#' @param bh_boundary Boundary window type for BH diagnostics ("mult" or "add").
#' @param bh_win Relative window size for BH boundary support when bh_boundary="mult".
#' @param bh_delta Additive window size when bh_boundary="add".
#'
#' @return A list with elements \code{tables}, \code{plots}, \code{objects}, and \code{params}.
#' @export
dfdr_run_all <- function(xs,
                         alpha_main = 0.01,
                         alphas = c(1e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1),
                         eps = 0.2,
                         win_rel = 0.2,
                         truncation = "warn_drop",
                         k_fixed = 10,
                         k_sqrt_mult = 2,
                         qband_breaks = c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5),
                         low_conf = c(0.2, 0.5),
                         min_N_equalchance = 2000,
                         pep_binwidth = 0.05,
                         pep_n_min = 200,
                         pep_max = 0.5,
                         # ---- NEW: PEP sanity + alt reliability (TDC-style)
                         pep_sanity_thresholds = c(0.01, 0.05, 0.1, 0.2),
                         pep_density_max = 0.5,
                         pep_rel_tdc_breaks = seq(0, 0.5, by = 0.05),
                         pep_rel_tdc_add_decoy = 0L,
                         # ---- p-value diagnostics defaults
                         compute_pseudo_pvalues = TRUE,
                         pseudo_pvalue_method = "decoy_ecdf",
                         p_u_grid = c(seq(0.001, 0.1, by = 0.001), seq(0.11, 1, by = 0.01)),
                         p_u_max = 0.1,
                         p_stratify = NULL,
                         p_min_n_calib = 200,
                         p_lambdas = seq(0.5, 0.95, by = 0.05),
                         p_min_n_pi0 = 2000,
                         bh_boundary = c("mult", "add"),
                         bh_win = 0.2,
                         bh_delta = NULL) {

  bh_boundary <- match.arg(bh_boundary)

  if (!is.list(xs) || is.null(names(xs)) || any(names(xs) == "")) {
    rlang::abort("`xs` must be a *named* list of dfdr_tbl objects.")
  }
  purrr::walk(xs, validate_dfdr_tbl)

  # ---- Add pseudo p-values if requested
  if (compute_pseudo_pvalues) {
    xs <- add_pseudo_pvalues_to_list(xs, method = pseudo_pvalue_method)
  }

  # ---- Core target/decoy tables (stacked with list label)
  headline <- purrr::imap_dfr(xs, ~dplyr::mutate(dfdr_headline(.x, alpha_main, k_fixed, k_sqrt_mult), list = .y))

  # ----

  stability <- purrr::imap_dfr(xs, ~dplyr::mutate(dfdr_curve_stability(.x, alphas, k_fixed, k_sqrt_mult), list = .y))
  dwin <- purrr::imap_dfr(xs, ~dplyr::mutate(dfdr_local_tail(.x, alphas, win_rel, truncation), list = .y))
  elasticity <- purrr::imap_dfr(xs, ~dplyr::mutate(dfdr_elasticity(.x, alphas, eps), list = .y))

  # ---- Calibration: PEP reliability and sumpep (only for those with PEP)
  pep_rel <- list()
  sumpep <- list()

  # ---- PEP sanity + alternative (TDC-style) reliability
  pep_sanity <- list()
  pep_rel_tdc <- list()

  # ---- Equal-chance (q-band) diagnostics (always computed)
  eqchance <- list()

  # ---- p-value diagnostics (only if p exists and has finite values)
  p_calib <- list()
  p_pi0 <- list()
  bh_head <- list()
  bh_elast <- list()

  for (nm in names(xs)) {
    x <- xs[[nm]]

    # --- PEP-based diagnostics
    has_pep <- ("pep" %in% names(x)) && any(is.finite(x$pep))

    if (has_pep) {
      pr <- dfdr_pep_reliability(x, binwidth = pep_binwidth, n_min = pep_n_min, pep_max = pep_max)
      pep_rel[[nm]] <- pr
      sumpep[[nm]] <- dfdr_sumpep(x, alphas)

      # Sanity: decoys getting small PEP?
      pep_sanity[[nm]] <- dfdr_pep_decoy_sanity(x, thresholds = pep_sanity_thresholds)

      # Alternative reliability: D/T proxy within PEP bins
      pep_rel_tdc[[nm]] <- dfdr_pep_reliability_tdc(
        x,
        breaks = pep_rel_tdc_breaks,
        add_decoy = pep_rel_tdc_add_decoy
      )
    }

    # --- Equal-chance by q-bands
    ec <- dfdr_equal_chance_qbands(x, breaks = qband_breaks, low_conf = low_conf, min_N = min_N_equalchance)
    eqchance[[nm]] <- ec

    # --- p-value diagnostics (only if p exists and usable)
    if ("p" %in% names(x) && any(is.finite(x$p))) {
      p_calib[[nm]] <- dfdr_p_calibration(
        x,
        u_grid = p_u_grid,
        u_max = p_u_max,
        stratify = p_stratify,
        min_n = p_min_n_calib
      )

      p_pi0[[nm]] <- dfdr_pi0_storey(
        x,
        lambdas = p_lambdas,
        stratify = p_stratify,
        min_n = p_min_n_pi0,
        clamp = TRUE
      )

      bh_head[[nm]] <- dfdr_bh_diagnostics(
        x,
        alpha = alpha_main,
        boundary = bh_boundary,
        win = bh_win,
        delta = bh_delta
      )

      bh_elast[[nm]] <- dfdr_bh_elasticity(
        x,
        alphas = alphas,
        eps = eps
      )

    }
  }

  # ---- merge additional columns into headline before flagging

  # Merge D_alpha_win from local tail (at alpha_main)
  dwin_main <- dwin |>
    dplyr::filter(alpha == alpha_main) |>
    dplyr::select(list, D_alpha_win)

  headline <- headline |>
    dplyr::left_join(dwin_main, by = "list")

  # Merge IPE from PEP reliability
  if (length(pep_rel) > 0) {
    ipe_tbl <- purrr::imap_dfr(pep_rel, ~tibble::tibble(list = .y, IPE = .x$IPE))
    headline <- headline |>
      dplyr::left_join(ipe_tbl, by = "list")
  }

  # Merge effect_abs from equal-chance pooled
  if (length(eqchance) > 0) {
    effect_tbl <- purrr::imap_dfr(eqchance, ~dplyr::mutate(.x$pooled, list = .y)) |>
      dplyr::select(.data$list, .data$effect_abs)
    headline <- headline |> dplyr::left_join(effect_tbl, by = "list")
  }

  # ---- flag (with columns present but possibly NA)
  headline <- flag_headline(headline)

  # ---- Core plots (common plots for multi-table)
  plots <- list(
    dalpha = dfdr_plot_dalpha(stability),
    cv = dfdr_plot_cv(stability),
    dwin = dfdr_plot_dwin(dwin[!is.na(dwin$D_alpha_win), , drop = FALSE]),
    elasticity = dfdr_plot_elasticity(elasticity)
  )

  # ---- Scope disagreement matrix (if multiple lists provided)
  if (length(xs) > 1) {
    plots$scope_disagreement <- dfdr_plot_scope_disagreement_matrix(xs, alpha_main)
  }

  # ---- Score distribution plots per list
  for (nm in names(xs)) {
    plots[[paste0("score_dist__", nm)]] <-
      dfdr_plot_score_distributions(xs[[nm]], title = paste0("Score distributions: ", nm))
  }

  # ---- PEP reliability plots per list (if available)
  for (nm in names(pep_rel)) {
    plots[[paste0("pep_reliability__", nm)]] <-
      dfdr_plot_pep_reliability(pep_rel[[nm]]$bins, title = paste0("PEP reliability: ", nm))
  }

  # ---- PEP density sanity plots (targets vs decoys)
  for (nm in names(pep_sanity)) {
    plots[[paste0("pep_density__", nm)]] <-
      dfdr_plot_pep_density_by_decoy(
        xs[[nm]],
        pep_max = pep_density_max,
        title = paste0("PEP density (targets vs decoys): ", nm)
      )
  }

  # ---- PEP reliability (TDC-style) plots per list
  for (nm in names(pep_rel_tdc)) {
    plots[[paste0("pep_reliability_tdc__", nm)]] <-
      dfdr_plot_pep_reliability_tdc(
        pep_rel_tdc[[nm]],
        title = paste0("PEP reliability (TDC-style): ", nm)
      )
  }

  # ---- Equal-chance plots per list
  for (nm in names(eqchance)) {
    plots[[paste0("equal_chance__", nm)]] <-
      dfdr_plot_equal_chance(eqchance[[nm]]$bands, title = paste0("Equal-chance by q-band: ", nm))
  }

  # ---- p-value / pseudo-p-value density plots per list (if p exists)
  for (nm in names(xs)) {
    x <- xs[[nm]]
    if ("p" %in% names(x) && any(is.finite(x$p))) {
      plots[[paste0("p_density__", nm)]] <-
        dfdr_plot_p_density_by_decoy(
          x,
          p_max = 1,
          log10_x = FALSE,
          title = paste0("(Pseudo-)p-value density (targets vs decoys): ", nm)
        )
    }
  }

  # ---- p-value plots per list (if available)
  for (nm in names(p_calib)) {
    x<-xs[[nm]]
    if ("p" %in% names(x) && any(is.finite(x$p))) {
    plots[[paste0("Calibration_pvalues__", nm)]] <-
      dfdr_plot_p_calibration2(x,sel="all")
    }
  }
  for (nm in names(p_calib)) {
    x<-xs[[nm]]
    if ("p" %in% names(x) && any(is.finite(x$p))) {
      plots[[paste0("Calibration_pvalues_decoys__", nm)]] <-
        dfdr_plot_p_calibration2(x,sel="decoy")
    }
  }
  for (nm in names(p_calib)) {
    x<-xs[[nm]]
    if ("p" %in% names(x) && any(is.finite(x$p))) {
      plots[[paste0("Calibration_pvalues_targets__", nm)]] <-
        dfdr_plot_p_calibration2(x,sel="target")
    }
  }
  for (nm in names(p_calib)) {
    x<-xs[[nm]]
    if ("p" %in% names(x) && any(is.finite(x$p))) {
      plots[[paste0("Estimated_fdr__", nm)]] <-
        dfdr_plot_fdrhat_pi0(x,sel="all")$plot
    }
  }
  for (nm in names(p_calib)) {
    x<-xs[[nm]]
    if ("p" %in% names(x) && any(is.finite(x$p))) {
      plots[[paste0("Estimated_fdr_decoys__", nm)]] <-
        dfdr_plot_fdrhat_pi0(x,sel="decoy")$plot
    }
  }
  for (nm in names(p_calib)) {
    x<-xs[[nm]]
    if ("p" %in% names(x) && any(is.finite(x$p))) {
      plots[[paste0("Estimated_fdr_targets__", nm)]] <-
        dfdr_plot_fdrhat_pi0(x,sel="target")$plot
    }
  }
  # for (nm in names(p_calib)) {
  #   plots[[paste0("Super_uniformity_pvalues__", nm)]] <-
  #     dfdr_plot_p_calibration(p_calib[[nm]]$ecdf, title = paste0("P-value calibration: ", nm))
  # }
  for (nm in names(p_pi0)) {
    plots[[paste0("pi0_Storey__", nm)]] <-
      dfdr_plot_pi0(p_pi0[[nm]]$pi0, title = paste0("Storey pi0(lambda): ", nm))
  }
  for (nm in names(bh_elast)) {
    plots[[paste0("BH_elasticity__", nm)]] <-
      dfdr_plot_bh_elasticity(bh_elast[[nm]], title = paste0("BH elasticity: ", nm))
  }

  # ---- Stack BH headline tables (if any)
  bh_headline_tbl <- NULL
  if (length(bh_head)) {
    bh_headline_tbl <- purrr::imap_dfr(bh_head, ~dplyr::mutate(.x$headline, list = .y))
  } else {
    bh_headline_tbl <- tibble::tibble()
  }

  # ---- Stack BH elasticity tables (if any)
  bh_elasticity_tbl <- NULL
  if (length(bh_elast)) {
    bh_elasticity_tbl <- purrr::imap_dfr(bh_elast, ~dplyr::mutate(.x, list = .y))
  } else {
    bh_elasticity_tbl <- tibble::tibble()
  }

  list(
    tables = list(
      # target/decoy
      headline = headline,  # contains flags
      stability = stability,
      local_tail = dwin,
      elasticity = elasticity,

      # PEP
      sumpep = sumpep,
      pep_reliability_bins = purrr::imap(pep_rel, ~.x$bins),
      pep_IPE = purrr::imap_dfr(pep_rel, ~tibble::tibble(list = .y, IPE = .x$IPE)),

      # PEP sanity + alt reliability
      pep_decoy_sanity = purrr::imap_dfr(pep_sanity, ~dplyr::mutate(.x, list = .y)),
      pep_reliability_tdc = purrr::imap_dfr(pep_rel_tdc, ~dplyr::mutate(.x, list = .y)),

      # equal-chance
      equal_chance_bands = purrr::imap(eqchance, ~.x$bands),
      equal_chance_pooled = purrr::imap_dfr(eqchance, ~dplyr::mutate(.x$pooled, list = .y)),

      # p-value calibration + pi0
      p_calibration_ecdf = purrr::imap(p_calib, ~.x$ecdf),
      p_calibration_summary = purrr::imap_dfr(p_calib, ~dplyr::mutate(.x$summary, list = .y)),
      pi0_storey = purrr::imap(p_pi0, ~.x$pi0),
      pi0_storey_summary = purrr::imap_dfr(p_pi0, ~dplyr::mutate(.x$summary, list = .y)),

      # BH
      bh_headline = bh_headline_tbl,
      bh_elasticity = bh_elasticity_tbl
    ),
    plots = plots,
    objects = xs,
    params = list(
      alpha_main = alpha_main,
      alphas = alphas,
      eps = eps,
      win_rel = win_rel,
      truncation = truncation,
      # PEP params
      pep_sanity_thresholds = pep_sanity_thresholds,
      pep_density_max = pep_density_max,
      pep_rel_tdc_breaks = pep_rel_tdc_breaks,
      pep_rel_tdc_add_decoy = pep_rel_tdc_add_decoy,
      # p-value params
      p_u_max = p_u_max,
      p_stratify = p_stratify,
      p_min_n_calib = p_min_n_calib,
      p_lambdas = p_lambdas,
      p_min_n_pi0 = p_min_n_pi0,
      bh_boundary = bh_boundary,
      bh_win = bh_win,
      bh_delta = bh_delta
    )
  )
}

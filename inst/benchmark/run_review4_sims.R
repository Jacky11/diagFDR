# Additional simulations of the accompanying article.
#   part "grid"      : the 540 datasets of the main grid (regenerated exactly as in run_simulation_v2.R) with the
#                      per-site paired single model and the consensus single model restricted to complete pairs
#   part "n2"        : two samples per condition (rho 0/0.5/0.9, three MCAR shares, ten seeds)
#   part "scenarios" : departures from the main simulation model (site clustering on proteins, condition-specific
#                      correlation, heterogeneous variances, unbalanced design and missingness, unmodelled batches,
#                      PTM loading offsets with normalisation, outlying values)
#   part "hetero"    : variances heterogeneous over sites (moderation prior df of 2-3, as in the real datasets)
#   part "variance"  : site residual variance several times the protein's, as in the real datasets
#   part "negative"  : negative within-sample correlation
#   cfg$seeds        : seeds per setting for the parts hetero, variance and negative (1:5 by default; 1:20 in the article)
# Usage, from the package folder: devtools::load_all(); review4_config <- list(part = "n2", out = "..."); source("inst/benchmark/run_review4_sims.R")
local({
  cfg <- utils::modifyList(list(part = "n2", out = file.path(normalizePath("."), "benchmark", "results", "review4"), cores = 2, seeds = 1:5,
                                pkg_dir = normalizePath(".")),
                           if (exists("review4_config", envir = globalenv())) get("review4_config", envir = globalenv()) else list())
  out <- file.path(cfg$out, cfg$part); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (cfg$part == "grid") {
    grid <- expand.grid(seed = 1:10, rho = c(0, 0.5, 0.9), pi_mcar = c(1, 0.5, 0.2), rate = c(0.1, 0.2, 0.4), n_rep = c(3, 5), KEEP.OUT.ATTRS = FALSE)
    grid$scenario <- "grid"
  } else if (cfg$part == "n2") {
    grid <- expand.grid(seed = 1:10, rho = c(0, 0.5, 0.9), pi_mcar = c(1, 0.5, 0.2), rate = 0.2, n_rep = 2, KEEP.OUT.ATTRS = FALSE)
    grid$scenario <- "n2"
  } else if (cfg$part == "variance") {
    # site residual variance several times the protein's, as in the real datasets (scripts/review4_variance_ratio.R),
    # alone and with heterogeneous variances over sites
    grid <- rbind(expand.grid(seed = cfg$seeds, rho = c(0, 0.5, 0.9), pi_mcar = 0.5, rate = 0.2, n_rep = c(2, 3, 5), scenario = "site_var", KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE),
                  expand.grid(seed = cfg$seeds, rho = c(0, 0.5, 0.9), pi_mcar = 0.5, rate = 0.2, n_rep = c(3, 5), scenario = "site_var_hetero", KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  } else if (cfg$part == "negative") {
    # negative within-sample correlation (enriched and total fractions not from exactly the same material)
    grid <- expand.grid(seed = cfg$seeds, rho = c(-0.3, -0.6), pi_mcar = 0.5, rate = 0.2, n_rep = c(3, 5), KEEP.OUT.ATTRS = FALSE)
    grid$scenario <- "negative"
  } else if (cfg$part == "hetero") {
    # realistic variance heterogeneity: sd_heterogeneity = 1 gives moderated degrees of freedom close to those of the real
    # datasets (prior df of about 2-3), whereas the default simulation has an almost infinite prior
    grid <- expand.grid(seed = cfg$seeds, rho = c(0, 0.5, 0.9), pi_mcar = 0.5, rate = 0.2, n_rep = c(2, 3, 5), KEEP.OUT.ATTRS = FALSE)
    grid$scenario <- "hetero"
  } else {
    sc <- c("baseline", "clustering", "rho_by_condition", "heterogeneous_variance", "unbalanced_missingness", "batch",
            "loading_normalised", "outliers_2pc", "outliers_5pc")
    sc <- c(sc, "loading_global", "loading_none")
    if (!is.null(cfg$scenarios)) sc <- intersect(sc, cfg$scenarios)
    grid <- rbind(expand.grid(seed = 1:5, rho = 0.5, pi_mcar = 0.5, rate = 0.2, n_rep = c(3, 5), scenario = sc, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE),
                  expand.grid(seed = 1:10, rho = 0.5, pi_mcar = 0.5, rate = 0.2, n_rep = 35, scenario = intersect("unbalanced_3v5", c(sc, if (is.null(cfg$scenarios)) "unbalanced_3v5")), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  }
  if (!is.null(cfg$only)) grid <- grid[cfg$only, ]
  grid$id <- sprintf("%s_rep%d_rate%02d_mcar%03d_rho%02d_seed%02d", grid$scenario, grid$n_rep, round(100 * grid$rate), round(100 * grid$pi_mcar), round(10 * grid$rho), grid$seed)
  run_one <- function(g, out, part) {
    f <- file.path(out, paste0(g$id, ".rds")); if (file.exists(f)) return(f)
    opt <- switch(g$scenario,
      clustering = list(site_clustering = 0.06), rho_by_condition = list(rho_condition = c(0.9, 0.1)),
      heterogeneous_variance = list(sd_heterogeneity = 0.5), unbalanced_missingness = list(rate_B_factor = 2),
      batch = list(batch_sd = c(0.15, 0.3)), loading_normalised = list(loading_sd = 0.3), loading_global = list(loading_sd = 0.3), loading_none = list(loading_sd = 0.3),
      outliers_2pc = list(outlier_frac = 0.02), outliers_5pc = list(outlier_frac = 0.05), hetero = list(sd_heterogeneity = 1),
      site_var = list(site_sd_factor = 2.5), site_var_hetero = list(site_sd_factor = 2.5, sd_heterogeneity = 0.5), list())
    n_rep <- if (g$n_rep == 35) c(3, 5) else g$n_rep
    set.seed(g$seed)
    sim <- do.call(simulate_ptm_features, c(list(n_rep = n_rep, rho = g$rho, rate = c(g$rate / 2, g$rate), pi_mcar = g$pi_mcar, seed = g$seed), opt))
    x0 <- filter_min_observed(sim_to_experiment(sim), min_n = 2)
    if (g$scenario == "loading_normalised") x0 <- normalise(x0, c("protein", "ptm"), "median")
    if (g$scenario == "loading_global") x0 <- normalise(x0, c("protein", "ptm"), "global_median")
    truth <- truth_table(x0); ctr <- make_contrasts(x0); x <- match_ptm_to_protein(x0, "id")
    xc <- x; na <- is.na(x$ptm$intensity) | is.na(x$ptm$protein_intensity); xc$ptm$intensity[na] <- NA; xc$ptm$protein_intensity[na] <- NA
    q <- function(e) suppressWarnings(suppressMessages(e))
    L <- function(xx, ...) q(method_limma(xx, ctr, analyses = "ptm_DPU", pi0_method = "bh", ...))
    arms <- if (part == "grid") list(per_site = function() L(x, blocking = "fixed"), consensus_complete = function() L(xc, blocking = "random")) else c(list(
      consensus = function() L(x, blocking = "random", complete_pairs = FALSE), consensus_complete = function() L(xc, blocking = "random"),
      per_site = function() L(x, blocking = "fixed"), unpaired = function() L(x, pair_by_sample = FALSE),
      unpaired_ordinary = function() L(x, pair_by_sample = FALSE, moderate = FALSE),
      msqrob2ptm = function() q(method_msqrob2ptm(x, ctr, analyses = "ptm_DPU")),
      no_correction = function() { r <- q(method_limma(x, ctr, analyses = "ptm_DPA", pi0_method = "bh")); r$analysis <- "ptm_DPU"; r }),
      if (g$scenario %in% c("baseline", "outliers_2pc", "outliers_5pc")) list(msqrob2ptm_nonrobust = function() q(method_msqrob2ptm(x, ctr, analyses = "ptm_DPU", robust = FALSE)))
      else if ((part == "n2" && g$pi_mcar == 0.5 && g$seed <= 3) || (part == "hetero" && g$seed <= 3)) list(msstatsptm = function() q(method_msstatsptm(x, ctr, analyses = "ptm_DPU"))))
    arms <- Filter(Negate(is.null), arms); res <- list(); rt <- list()
    for (a in names(arms)) {
      t0 <- Sys.time(); r <- tryCatch(arms[[a]](), error = function(e) NULL); rt[[a]] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (!is.null(r) && nrow(r)) { r$method <- a; res[[a]] <- r[r$analysis == "ptm_DPU", ] }
    }
    res <- do.call(rbind, res)
    meta <- data.frame(id = g$id, scenario = g$scenario, seed = g$seed, n_rep = g$n_rep, rate = g$rate, pi_mcar = g$pi_mcar, rho = g$rho)
    cls <- table(sim$truth$class)
    out_l <- list(metrics = cbind(meta, evaluate_against_truth(res, truth, c(0.01, 0.05, 0.1))),
                  by_class = cbind(meta, evaluate_by_class(res, truth, c(0.01, 0.05, 0.1))),
                  runtime = cbind(meta, data.frame(arm = names(rt), seconds = unlist(rt))),
                  classes = cbind(meta, data.frame(class = names(cls), n = as.integer(cls))),
                  fdp_curve = cbind(meta, fdp_curve(res, truth)),
                  df = cbind(meta, stats::aggregate(df ~ method, data = res[is.finite(res$df), ], FUN = stats::median)))
    if (part == "grid" && g$rate == 0.2 && g$pi_mcar == 0.5 && g$seed <= 5) out_l$sites <- res   # the 30 datasets of the combination analysis
    saveRDS(out_l, f); f
  }
  rows <- split(grid, seq_len(nrow(grid)))
  if (cfg$cores > 1) {
    cl <- parallel::makeCluster(cfg$cores)
    parallel::clusterCall(cl, function(d) { suppressMessages(devtools::load_all(d, quiet = TRUE)); data.table::setDTthreads(1) }, cfg$pkg_dir)
    parallel::parLapplyLB(cl, rows, run_one, out = out, part = cfg$part); parallel::stopCluster(cl)
  } else lapply(rows, run_one, out = out, part = cfg$part)
  all <- lapply(file.path(out, paste0(grid$id, ".rds")), readRDS)
  for (p in c("metrics", "by_class", "runtime", "classes", "fdp_curve", "df"))
    data.table::fwrite(data.table::rbindlist(lapply(all, `[[`, p), fill = TRUE), file.path(cfg$out, paste0("review4_", cfg$part, "_", p, ".tsv")), sep = "\t")
  message("done: ", nrow(grid), " datasets (", cfg$part, ")")
})

# The consensus single model restricted to complete pairs on the datasets of run_sensitivity.R
# (unit coefficient, heterogeneous correlations, summarisation filter), regenerated identically.
# Usage, from the package folder: devtools::load_all(); sens_cc_config <- list(out = "..."); source("inst/benchmark/run_sensitivity_complete_pairs.R")
local({
  cfg <- utils::modifyList(list(out = file.path(normalizePath("."), "benchmark", "results", "sensitivity_cc"), seeds = 1:5, cores = 2, pkg_dir = normalizePath(".")),
                           if (exists("sens_cc_config", envir = globalenv())) get("sens_cc_config", envir = globalenv()) else list())
  dir.create(cfg$out, recursive = TRUE, showWarnings = FALSE)
  A <- expand.grid(seed = cfg$seeds, beta = c(0.6, 0.8, 1, 1.2, 1.5), n_rep = c(3, 5), rho = 0.5, pi_mcar = 0.5, rate = 0.2, het = FALSE, filter = TRUE, part = "A", stringsAsFactors = FALSE)
  B <- expand.grid(seed = cfg$seeds, beta = 1, n_rep = c(3, 5), rho = 0.45, pi_mcar = c(1, 0.2), rate = 0.2, het = c(FALSE, TRUE), filter = TRUE, part = "B", stringsAsFactors = FALSE)
  C <- expand.grid(seed = cfg$seeds, beta = 1, n_rep = c(3, 5), rho = 0.9, pi_mcar = c(1, 0.2), rate = 0.4, het = FALSE, filter = c(TRUE, FALSE), part = "C", stringsAsFactors = FALSE)
  grid <- rbind(A, B, C)
  grid$id <- sprintf("%s_beta%03d_rep%d_rho%02d_mcar%03d_rate%02d_het%d_filt%d_seed%02d", grid$part, round(100 * grid$beta), grid$n_rep,
                     round(100 * grid$rho), round(100 * grid$pi_mcar), round(100 * grid$rate), grid$het, grid$filter, grid$seed)
  run_one <- function(g, out) {
    f <- file.path(out, paste0(g$id, ".rds")); if (file.exists(f)) return(f)
    set.seed(g$seed)                                          # as in run_sensitivity.R: the per-site rho is drawn after this
    rho <- if (g$het) stats::runif(1500, 0, 0.9) else g$rho
    sim <- simulate_ptm_features(n_rep = g$n_rep, rho = rho, rate = c(g$rate / 2, g$rate), pi_mcar = g$pi_mcar, seed = g$seed, beta = g$beta)
    x0 <- filter_min_observed(sim_to_experiment(sim, condition_filter = g$filter), min_n = 2)
    truth <- truth_table(x0); ctr <- make_contrasts(x0); x <- match_ptm_to_protein(x0, "id")
    r <- suppressWarnings(method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh", blocking = "random", complete_pairs = TRUE)); r$method <- "paired_random_complete"
    r0 <- suppressWarnings(method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh", blocking = "random", complete_pairs = FALSE)); r0$method <- "paired_random_check"
    r <- rbind(r, r0)                                         # the check arm must reproduce paired_random of run_sensitivity.R
    meta <- data.frame(id = g$id, part = g$part, seed = g$seed, beta = g$beta, n_rep = g$n_rep, rho = g$rho, het = g$het, pi_mcar = g$pi_mcar, rate = g$rate, filter = g$filter)
    saveRDS(list(metrics = cbind(meta, evaluate_against_truth(r, truth, c(0.01, 0.05, 0.1))), by_class = cbind(meta, evaluate_by_class(r, truth, c(0.01, 0.05, 0.1)))), f); f
  }
  rows <- split(grid, seq_len(nrow(grid)))
  if (cfg$cores > 1) {
    cl <- parallel::makeCluster(cfg$cores)
    parallel::clusterCall(cl, function(d) { suppressMessages(devtools::load_all(d, quiet = TRUE)); data.table::setDTthreads(1) }, cfg$pkg_dir)
    parallel::parLapplyLB(cl, rows, run_one, out = cfg$out); parallel::stopCluster(cl)
  } else lapply(rows, run_one, out = cfg$out)
  all <- lapply(file.path(cfg$out, paste0(grid$id, ".rds")), readRDS)
  for (p in c("metrics", "by_class")) data.table::fwrite(data.table::rbindlist(lapply(all, `[[`, p)), file.path(dirname(cfg$out), paste0("sensitivity_cc_", p, ".tsv")), sep = "\t")
  message("done: ", nrow(grid), " datasets")
})

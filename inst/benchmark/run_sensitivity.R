# Benchmark 1c: sensitivity analyses.
#   A. unit-coefficient assumption: sites follow their protein with coefficient beta != 1
#   B. heterogeneous within-sample correlation: rho drawn per site from U(0, 0.9) vs constant 0.45
#   C. summarisation filter: with and without excluding features not observed in every condition
# Arms: single model paired with a fixed sample effect, paired with a random sample effect (consensus
# correlation), unpaired; msqrob2PTM; no correction. No imputation.
# Usage, from the package folder: devtools::load_all(); source("inst/benchmark/run_sensitivity.R")
local({
  cfg <- utils::modifyList(list(out = file.path(normalizePath("."), "benchmark", "results", "sensitivity"), seeds = 1:5, cores = 2,
                                pkg_dir = normalizePath(".")),
                           if (exists("sensitivity_config", envir = globalenv())) get("sensitivity_config", envir = globalenv()) else list())
  dir.create(cfg$out, recursive = TRUE, showWarnings = FALSE)
  A <- expand.grid(seed = cfg$seeds, beta = c(0.6, 0.8, 1, 1.2, 1.5), n_rep = c(3, 5), rho = 0.5, pi_mcar = 0.5, rate = 0.2, het = FALSE, filter = TRUE, part = "A", stringsAsFactors = FALSE)
  B <- expand.grid(seed = cfg$seeds, beta = 1, n_rep = c(3, 5), rho = 0.45, pi_mcar = c(1, 0.2), rate = 0.2, het = c(FALSE, TRUE), filter = TRUE, part = "B", stringsAsFactors = FALSE)
  C <- expand.grid(seed = cfg$seeds, beta = 1, n_rep = c(3, 5), rho = 0.9, pi_mcar = c(1, 0.2), rate = 0.4, het = FALSE, filter = c(TRUE, FALSE), part = "C", stringsAsFactors = FALSE)
  grid <- rbind(A, B, C)
  grid$id <- sprintf("%s_beta%03d_rep%d_rho%02d_mcar%03d_rate%02d_het%d_filt%d_seed%02d", grid$part, round(100 * grid$beta), grid$n_rep,
                     round(100 * grid$rho), round(100 * grid$pi_mcar), round(100 * grid$rate), grid$het, grid$filter, grid$seed)
  run_one <- function(g, out) {
    f <- file.path(out, paste0(g$id, ".rds")); if (file.exists(f)) return(f)
    set.seed(g$seed)
    rho <- if (g$het) stats::runif(1500, 0, 0.9) else g$rho
    sim <- simulate_ptm_features(n_rep = g$n_rep, rho = rho, rate = c(g$rate / 2, g$rate), pi_mcar = g$pi_mcar, seed = g$seed, beta = g$beta)
    x0 <- filter_min_observed(sim_to_experiment(sim, condition_filter = g$filter), min_n = 2)
    truth <- truth_table(x0); ctr <- make_contrasts(x0); x <- match_ptm_to_protein(x0, "id")
    n_feat <- c(protein = mean(x$protein$n_features, na.rm = TRUE), ptm = mean(x$ptm$n_features, na.rm = TRUE))
    res <- list()
    add <- function(label, r) { r$method <- label; res[[label]] <<- r }
    add("paired_random", suppressWarnings(method_limma(x, ctr, pi0_method = "bh")))
    add("paired_fixed", suppressWarnings(method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh", blocking = "fixed")))
    add("unpaired", suppressWarnings(method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh", pair_by_sample = FALSE)))
    nc <- res$paired_random[res$paired_random$analysis == "ptm_DPA", ]; nc$analysis <- "ptm_DPU"; add("no_correction", nc)
    if ("msqrob2ptm" %in% available_methods()) add("msqrob2ptm", suppressWarnings(suppressMessages(method_msqrob2ptm(x, ctr))))
    res <- do.call(rbind, res)
    meta <- data.frame(id = g$id, part = g$part, seed = g$seed, beta = g$beta, n_rep = g$n_rep, rho = g$rho, het = g$het, pi_mcar = g$pi_mcar,
                       rate = g$rate, filter = g$filter, mean_features_protein = n_feat[["protein"]], mean_features_ptm = n_feat[["ptm"]],
                       n_sites = nrow(x$ptm$intensity), n_proteins = nrow(x$protein$intensity))
    out_l <- list(metrics = cbind(meta, evaluate_against_truth(res, truth, c(0.01, 0.05, 0.1))),
                  by_class = cbind(meta, evaluate_by_class(res, truth, c(0.01, 0.05, 0.1))))
    saveRDS(out_l, f); f
  }
  rows <- split(grid, seq_len(nrow(grid)))
  if (cfg$cores > 1) {
    cl <- parallel::makeCluster(cfg$cores)
    parallel::clusterCall(cl, function(d) { suppressMessages(devtools::load_all(d, quiet = TRUE)); data.table::setDTthreads(1) }, cfg$pkg_dir)
    parallel::parLapplyLB(cl, rows, run_one, out = cfg$out); parallel::stopCluster(cl)
  } else lapply(rows, run_one, out = cfg$out)
  all <- lapply(file.path(cfg$out, paste0(grid$id, ".rds")), readRDS)
  for (p in c("metrics", "by_class"))
    utils::write.table(do.call(rbind, lapply(all, `[[`, p)), file.path(dirname(cfg$out), paste0("sensitivity_", p, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  message("done: ", nrow(grid), " datasets")
})

# Benchmark 1d: combining the strategies' calls (intersection, union, majority) — does it improve
# error control or sensitivity over the single strategies? Per-site results of every strategy are
# kept so that any combination rule can be evaluated against the truth.
# Grid (subset of the main grid): n_rep {3, 5} x 20 % missing x MCAR 0.5 x rho {0, 0.5, 0.9} x seeds 1-5.
# Usage, from the package folder: devtools::load_all(); source("inst/benchmark/run_combination.R")
local({
  cfg <- utils::modifyList(list(out = file.path(normalizePath("."), "benchmark", "results", "combination"), seeds = 1:5, cores = 2,
                                pkg_dir = normalizePath("."), rhos = c(0, 0.5, 0.9), n_rep = c(3, 5), rate = 0.2, pi_mcar = 0.5),
                           if (exists("combination_config", envir = globalenv())) get("combination_config", envir = globalenv()) else list())
  if (!nzchar(cfg$pkg_dir)) cfg$pkg_dir <- normalizePath(".")
  dir.create(cfg$out, recursive = TRUE, showWarnings = FALSE)
  grid <- expand.grid(seed = cfg$seeds, rho = cfg$rhos, n_rep = cfg$n_rep, rate = cfg$rate, pi_mcar = cfg$pi_mcar, KEEP.OUT.ATTRS = FALSE)
  grid$id <- sprintf("rep%d_rate%02d_mcar%03d_rho%02d_seed%02d", grid$n_rep, round(100 * grid$rate), round(100 * grid$pi_mcar), round(10 * grid$rho), grid$seed)
  run_one <- function(g, out) {
    f <- file.path(out, paste0(g$id, ".rds")); if (file.exists(f)) return(f)
    set.seed(g$seed)
    sim <- simulate_ptm_features(n_rep = g$n_rep, rho = g$rho, rate = c(g$rate / 2, g$rate), pi_mcar = g$pi_mcar, seed = g$seed)
    x0 <- filter_min_observed(sim_to_experiment(sim), min_n = 2); truth <- truth_table(x0); ctr <- make_contrasts(x0)
    x <- match_ptm_to_protein(x0, "id")
    res <- list()
    add <- function(label, r) { r$method <- label; res[[label]] <<- r[r$analysis == "ptm_DPU", ] }
    add("limma", suppressWarnings(method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh")))
    add("limma_fixed", suppressWarnings(method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh", blocking = "fixed")))
    add("limma_unpaired", suppressWarnings(method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh", pair_by_sample = FALSE)))
    if ("msqrob2ptm" %in% available_methods()) add("msqrob2ptm", suppressWarnings(suppressMessages(method_msqrob2ptm(x, ctr, analyses = "ptm_DPU"))))
    if ("msstatsptm" %in% available_methods()) add("msstatsptm", suppressWarnings(suppressMessages(method_msstatsptm(x, ctr))))
    res <- do.call(rbind, res); rownames(res) <- NULL
    meta <- data.frame(id = g$id, seed = g$seed, n_rep = g$n_rep, rate = g$rate, pi_mcar = g$pi_mcar, rho = g$rho)
    saveRDS(list(meta = meta, res = res, truth = truth), f); f
  }
  rows <- split(grid, seq_len(nrow(grid)))
  if (cfg$cores > 1) {
    cl <- parallel::makeCluster(cfg$cores)
    parallel::clusterCall(cl, function(d) { suppressMessages(devtools::load_all(d, quiet = TRUE)); data.table::setDTthreads(1) }, cfg$pkg_dir)
    parallel::parLapplyLB(cl, rows, run_one, out = cfg$out); parallel::stopCluster(cl)
  } else lapply(rows, run_one, out = cfg$out)
  message("done: ", nrow(grid), " datasets in ", cfg$out)
})

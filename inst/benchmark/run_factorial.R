# Benchmark 1b: pairing x moderation factorial with the same estimator (limma single model), and
# fixed-effect vs consensus-correlation (random-effect) sample term, on regenerated datasets of the
# main simulation grid (no imputation). Each arm is a variant of method_limma():
#   paired_fixed      sample as a fixed effect (exact within-sample pairing, Eq. 6)
#   paired_random     sample as a random effect, consensus correlation (limma::duplicateCorrelation)
#   unpaired          no sample term (Eq. 5)
# each moderated (eBayes, robust) or not (ordinary t on the feature's own residual variance).
# Usage, from the package folder: devtools::load_all(); source("inst/benchmark/run_factorial.R")
# Settings: factorial_config <- list(out = "...", seeds = 1:10, rates = c(0.2), cores = 4)
local({
  cfg <- utils::modifyList(list(out = file.path(normalizePath("."), "benchmark", "results", "factorial"),
                                seeds = 1:10, rates = c(0.2, 0.4), pi_mcar = c(1, 0.5, 0.2), rhos = c(0, 0.5, 0.9),
                                n_rep = c(3, 5), cores = 4, pkg_dir = normalizePath(".")),
                           if (exists("factorial_config", envir = globalenv())) get("factorial_config", envir = globalenv()) else list())
  dir.create(cfg$out, recursive = TRUE, showWarnings = FALSE)
  grid <- expand.grid(seed = cfg$seeds, rho = cfg$rhos, pi_mcar = cfg$pi_mcar, rate = cfg$rates, n_rep = cfg$n_rep, KEEP.OUT.ATTRS = FALSE)
  grid$id <- sprintf("rep%d_rate%02d_mcar%03d_rho%02d_seed%02d", grid$n_rep, round(100 * grid$rate), round(100 * grid$pi_mcar), round(10 * grid$rho), grid$seed)
  run_one <- function(g, out) {
    f <- file.path(out, paste0(g$id, ".rds")); if (file.exists(f)) return(f)
    set.seed(g$seed)
    sim <- simulate_ptm_features(n_rep = g$n_rep, rho = g$rho, rate = c(g$rate / 2, g$rate), pi_mcar = g$pi_mcar, seed = g$seed)
    x0 <- filter_min_observed(sim_to_experiment(sim), min_n = 2); truth <- truth_table(x0); ctr <- make_contrasts(x0)
    x <- match_ptm_to_protein(x0, "id")
    arms <- list(
      paired_fixed_moderated    = list(blocking = "fixed",  moderate = TRUE,  pair_by_sample = TRUE),
      paired_fixed_ordinary     = list(blocking = "fixed",  moderate = FALSE, pair_by_sample = TRUE),
      paired_random_moderated   = list(blocking = "random", moderate = TRUE,  pair_by_sample = TRUE),
      paired_random_ordinary    = list(blocking = "random", moderate = FALSE, pair_by_sample = TRUE),
      unpaired_moderated        = list(blocking = "auto",   moderate = TRUE,  pair_by_sample = FALSE),
      unpaired_ordinary         = list(blocking = "auto",   moderate = FALSE, pair_by_sample = FALSE))
    res <- list(); times <- list()
    for (a in names(arms)) {
      t0 <- Sys.time()
      r <- suppressWarnings(suppressMessages(do.call(method_limma, c(list(x, ctr, analyses = "ptm_DPU", pi0_method = "bh"), arms[[a]]))))
      r$method <- a; res[[a]] <- r; times[[a]] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    }
    res <- do.call(rbind, res)
    meta <- data.frame(id = g$id, seed = g$seed, n_rep = g$n_rep, rate = g$rate, pi_mcar = g$pi_mcar, rho = g$rho)
    out_l <- list(metrics = cbind(meta, evaluate_against_truth(res, truth, c(0.01, 0.05, 0.1))),
                  by_class = cbind(meta, evaluate_by_class(res, truth, c(0.01, 0.05, 0.1))),
                  runtime = cbind(meta, data.frame(arm = names(times), seconds = unlist(times))))
    saveRDS(out_l, f); f
  }
  rows <- split(grid, seq_len(nrow(grid)))
  if (cfg$cores > 1) {
    cl <- parallel::makeCluster(cfg$cores)
    parallel::clusterCall(cl, function(d) { suppressMessages(devtools::load_all(d, quiet = TRUE)); data.table::setDTthreads(1) }, cfg$pkg_dir)
    parallel::parLapplyLB(cl, rows, run_one, out = cfg$out); parallel::stopCluster(cl)
  } else lapply(rows, run_one, out = cfg$out)
  all <- lapply(file.path(cfg$out, paste0(grid$id, ".rds")), readRDS)
  for (p in c("metrics", "by_class", "runtime"))
    utils::write.table(do.call(rbind, lapply(all, `[[`, p)), file.path(dirname(cfg$out), paste0("factorial_", p, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  message("done: ", nrow(grid), " datasets")
})

# Benchmark 1e: (a) MSstatsPTM (harmonised) on the seven seeds (4-10) not run in the main
# grid, for the 20 % missing-value settings, so that the headline comparisons rest on 30 datasets per
# entry for every strategy; (b) msqrob2PTM with and without robust estimation on ten datasets.
# Output rows use the schema of run_simulation_v2.R and are merged into sim2_*.tsv by
# draft/scripts/review_merge_extra.R.
# Usage, from the package folder: devtools::load_all(); source("inst/benchmark/run_extra_seeds.R")
local({
  cfg <- utils::modifyList(list(out = file.path(normalizePath("."), "benchmark", "results", "extra"), cores = 2,
                                seeds = 4:10, rates = 0.2, pkg_dir = normalizePath("."), fdr_levels = c(0.01, 0.05, 0.1)),
                           if (exists("extra_config", envir = globalenv())) get("extra_config", envir = globalenv()) else list())
  dir.create(cfg$out, recursive = TRUE, showWarnings = FALSE)
  grid <- expand.grid(seed = cfg$seeds, rho = c(0, 0.5, 0.9), pi_mcar = c(1, 0.5, 0.2), rate = cfg$rates, n_rep = c(3, 5), KEEP.OUT.ATTRS = FALSE)
  grid$part <- "msstats"
  rob <- expand.grid(seed = 1:10, rho = 0.5, pi_mcar = 0.5, rate = 0.2, n_rep = 3, KEEP.OUT.ATTRS = FALSE); rob$part <- "robust"
  grid <- rbind(grid, rob)
  grid$id <- sprintf("%s_rep%d_rate%02d_mcar%03d_rho%02d_seed%02d", grid$part, grid$n_rep, round(100 * grid$rate), round(100 * grid$pi_mcar), round(10 * grid$rho), grid$seed)
  run_one <- function(g, out, fdr_levels) {
    f <- file.path(out, paste0(g$id, ".rds")); if (file.exists(f)) return(f)
    set.seed(g$seed)
    sim <- simulate_ptm_features(n_rep = g$n_rep, rho = g$rho, rate = c(g$rate / 2, g$rate), pi_mcar = g$pi_mcar, seed = g$seed)
    x0 <- filter_min_observed(sim_to_experiment(sim), min_n = 2); truth <- truth_table(x0); ctr <- make_contrasts(x0)
    x <- match_ptm_to_protein(x0, "id")
    q <- function(e) suppressWarnings(suppressMessages(e))
    arms <- list(); times <- list()
    if (g$part == "msstats") {
      t0 <- Sys.time(); r <- q(method_msstatsptm(x, ctr)); times[["msstatsptm none"]] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      r$method <- "msstatsptm"; r$imputation <- "none"; arms[["msstatsptm none"]] <- r
    } else {
      t0 <- Sys.time(); r <- q(method_msqrob2ptm(x, ctr)); times[["msqrob2ptm none"]] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      r$method <- "msqrob2ptm"; r$imputation <- "none"; arms[["msqrob2ptm none"]] <- r
      t0 <- Sys.time(); r <- q(method_msqrob2ptm(x, ctr, robust = FALSE)); times[["msqrob2ptm_nonrobust none"]] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      r$method <- "msqrob2ptm_nonrobust"; r$imputation <- "none"; arms[["msqrob2ptm_nonrobust none"]] <- r
    }
    res <- do.call(rbind, arms); res$arm <- paste(res$method, res$imputation, sep = "|")
    ev_res <- res; ev_res$method <- ev_res$arm
    split_arm <- function(d) { d$imputation <- sub("^.*\\|", "", d$method); d$method <- sub("\\|.*$", "", d$method); d }
    meta <- data.frame(id = sub("^[a-z]+_", "", g$id), seed = g$seed, n_rep = g$n_rep, rate = g$rate, pi_mcar = g$pi_mcar, rho = g$rho)
    out_l <- list(metrics = cbind(meta, split_arm(evaluate_against_truth(ev_res, truth, fdr_levels))),
                  by_class = cbind(meta, split_arm(evaluate_by_class(ev_res, truth, fdr_levels))),
                  fdp = cbind(meta, split_arm(fdp_curve(ev_res, truth))),
                  runtime = cbind(meta, data.frame(arm = names(times), seconds = unlist(times))))
    saveRDS(out_l, f); f
  }
  rows <- split(grid, seq_len(nrow(grid)))
  if (cfg$cores > 1) {
    cl <- parallel::makeCluster(cfg$cores)
    parallel::clusterCall(cl, function(d) { suppressMessages(devtools::load_all(d, quiet = TRUE)); data.table::setDTthreads(1) }, cfg$pkg_dir)
    parallel::parLapplyLB(cl, rows, run_one, out = cfg$out, fdr_levels = cfg$fdr_levels); parallel::stopCluster(cl)
  } else lapply(rows, run_one, out = cfg$out, fdr_levels = cfg$fdr_levels)
  all <- lapply(file.path(cfg$out, paste0(grid$id, ".rds")), readRDS)
  for (p in c("metrics", "by_class", "fdp", "runtime"))
    utils::write.table(do.call(rbind, lapply(all, `[[`, p)), file.path(dirname(cfg$out), paste0("extra_", p, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  message("done: ", nrow(grid), " datasets")
})

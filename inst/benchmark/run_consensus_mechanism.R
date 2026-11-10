# Why does the consensus-correlation pairing exceed the nominal FDR at high correlation and
# missingness? Regenerates the worst factorial settings with per-site results and covariates, and tests
# consensus pairing restricted to complete site-protein pairs (unpaired values set to NA before fitting).
# Arms: consensus (random sample term), consensus on complete pairs, per-site (fixed sample term), unpaired.
# Usage, from the package folder: devtools::load_all(); mechanism_config <- list(out = "..."); source("inst/benchmark/run_consensus_mechanism.R")
local({
  cfg <- utils::modifyList(list(out = file.path(normalizePath("."), "benchmark", "results", "mechanism"), seeds = 1:10,
                                cores = 2, pkg_dir = normalizePath(".")),
                           if (exists("mechanism_config", envir = globalenv())) get("mechanism_config", envir = globalenv()) else list())
  dir.create(cfg$out, recursive = TRUE, showWarnings = FALSE)
  grid <- rbind(expand.grid(seed = cfg$seeds, n_rep = c(3, 5), rate = 0.4, pi_mcar = c(1, 0.2), rho = 0.9),
                expand.grid(seed = cfg$seeds, n_rep = c(3, 5), rate = 0.2, pi_mcar = 0.5, rho = c(0, 0.9)))
  grid$id <- sprintf("rep%d_rate%02d_mcar%03d_rho%02d_seed%02d", grid$n_rep, round(100 * grid$rate), round(100 * grid$pi_mcar), round(10 * grid$rho), grid$seed)
  run_one <- function(g, out) {
    f <- file.path(out, paste0(g$id, ".rds")); if (file.exists(f)) return(f)
    sim <- simulate_ptm_features(n_rep = g$n_rep, rho = g$rho, rate = c(g$rate / 2, g$rate), pi_mcar = g$pi_mcar, seed = g$seed)
    x0 <- filter_min_observed(sim_to_experiment(sim), min_n = 2); truth <- truth_table(x0); ctr <- make_contrasts(x0)
    x <- match_ptm_to_protein(x0, "id")
    S <- x$ptm$intensity; P <- x$ptm$protein_intensity; grp <- as.character(groups_of(x, "ptm"))
    xc <- x; na <- is.na(S) | is.na(P); xc$ptm$intensity[na] <- NA; xc$ptm$protein_intensity[na] <- NA
    q <- function(e) suppressWarnings(suppressMessages(e))
    res <- list(
      consensus = q(method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh", blocking = "random")),
      consensus_complete = q(method_limma(xc, ctr, analyses = "ptm_DPU", pi0_method = "bh", blocking = "random")),
      per_site = q(method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh", blocking = "fixed")),
      unpaired = q(method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh", pair_by_sample = FALSE)))
    for (a in names(res)) res[[a]]$method <- a
    res <- do.call(rbind, res)
    cov <- data.frame(feature_id = rownames(S),
                      complete_pairs = rowSums(!is.na(S) & !is.na(P)), site_only = rowSums(!is.na(S) & is.na(P)),
                      protein_only = rowSums(is.na(S) & !is.na(P)), n_site = rowSums(!is.na(S)),
                      min_pairs_group = apply(!is.na(S) & !is.na(P), 1, function(v) min(tapply(v, grp, sum))),
                      ptm_features = rowMeans(x$ptm$n_features, na.rm = TRUE),
                      obs_cor = vapply(seq_len(nrow(S)), function(i) {
                        cs <- S[i, ] - stats::ave(S[i, ], grp, FUN = function(v) mean(v, na.rm = TRUE))
                        cp <- P[i, ] - stats::ave(P[i, ], grp, FUN = function(v) mean(v, na.rm = TRUE))
                        ok <- !is.na(cs) & !is.na(cp); if (sum(ok) >= 4) suppressWarnings(stats::cor(cs[ok], cp[ok])) else NA_real_ }, 0))
    meta <- data.frame(id = g$id, seed = g$seed, n_rep = g$n_rep, rate = g$rate, pi_mcar = g$pi_mcar, rho = g$rho)
    tr <- truth[truth$analysis == "ptm_DPU", c("feature_id", "positive", "class")]
    sites <- merge(merge(res[res$analysis == "ptm_DPU", c("feature_id", "method", "log2FC", "se", "df", "pvalue", "adj_pvalue")], tr, by = "feature_id"), cov, by = "feature_id", all.x = TRUE)
    out_l <- list(meta = meta, metrics = cbind(meta, evaluate_against_truth(res, truth, c(0.01, 0.05, 0.1))), sites = cbind(meta[rep(1, nrow(sites)), c("id")], sites))
    names(out_l$sites)[1] <- "id"
    saveRDS(out_l, f); f
  }
  rows <- split(grid, seq_len(nrow(grid)))
  if (cfg$cores > 1) {
    cl <- parallel::makeCluster(cfg$cores)
    parallel::clusterCall(cl, function(d) { suppressMessages(devtools::load_all(d, quiet = TRUE)); data.table::setDTthreads(1) }, cfg$pkg_dir)
    parallel::parLapplyLB(cl, rows, run_one, out = cfg$out); parallel::stopCluster(cl)
  } else lapply(rows, run_one, out = cfg$out)
  all <- lapply(file.path(cfg$out, paste0(grid$id, ".rds")), readRDS)
  utils::write.table(do.call(rbind, lapply(all, `[[`, "metrics")), file.path(dirname(cfg$out), "mechanism_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  data.table::fwrite(data.table::rbindlist(lapply(all, function(z) cbind(z$meta[rep(1, nrow(z$sites)), c("n_rep", "rate", "pi_mcar", "rho")], z$sites))),
                     file.path(dirname(cfg$out), "mechanism_sites.tsv.gz"), sep = "\t")
  message("done: ", nrow(grid), " datasets")
})

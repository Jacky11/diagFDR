# Benchmark 1 (v2): feature-level simulations with imp4p missing-value mechanism.
#
# Grid: replicates {3, 5} x PTM missing rate {10, 20, 40 %} (protein rate = half) x
#       MCAR fraction {1, 0.5, 0.2} x within-sample site-protein correlation rho {0, 0.5, 0.9}
#       x seeds 1..10  ->  540 simulated datasets.
# Arms per dataset (DiffPTM median-polish summaries of features observed in every condition, no
# normalisation since the simulation has no technical bias, features kept if >= 2 observations in a group):
#   imputation none / imp4p MLE / imp4p SLSA  x  limma paired single model, msqrob2PTM
#   imputation none only: limma unpaired single model (DiffPTM v1), limma with precision
#                         weights (number of features per summary; exploratory), MSstatsPTM harmonised,
#                         no correction (PTM abundance used as usage)
#   feature-level input:  MSstatsPTM native with AFT imputation, and without (MBimpute = FALSE)
# MSstatsPTM arms (about 1 min each) run on `msstats_seeds` only.
#
# Light on the computer:
#   * few workers by default (`cores`, default 4), each limited to one thread (data.table, BLAS)
#     and run at "below normal" priority, so the computer stays usable;
#   * datasets are processed in batches; workers are restarted after each batch, which returns
#     their memory to the system;
#   * every dataset is saved as soon as it is done (one .rds file), and the summary tables are
#     rewritten after each batch, so partial results are usable and an interrupted run (even a
#     shutdown) resumes where it stopped: just source the script again;
#   * progress, elapsed time and estimated time left are printed after each batch and appended to
#     `progress.log` in the output folder.
#
# Output (in benchmark/results/):
#   sim2/<dataset>.rds    one file per simulated dataset
#   sim2_metrics.tsv      one row per dataset x arm x analysis x nominal FDR
#   sim2_by_class.tsv     call rate per truth class
#   sim2_fdp_curve.tsv    observed FDP over nominal FDR 0.1-20 %
#   sim2_runtime.tsv      seconds per arm and dataset
#
# Usage, from the package folder (preferably in a fresh R session):
#   devtools::load_all(); source("inst/benchmark/run_simulation_v2.R")
# To change settings, define a list `sim2_config` first, for example
#   sim2_config <- list(cores = 2)                 # lighter still
#   sim2_config <- list(cores = 6, batch_size = 12) # faster, if the computer copes
#   sim2_config <- list(max_batches = 1)           # run one batch, then stop (resume later)
# Variables left in the session by other scripts (n_seeds, out_dir, ...) are ignored.

local({
  defaults <- list(
    pkg_dir       = normalizePath("."),
    out_dir       = NULL,                 # default: <pkg_dir>/benchmark/results/sim2
    import_from   = NULL,                 # default: <pkg_dir>/benchmark/results (files of an earlier run)
    n_seeds       = 10,
    msstats_seeds = 1:3,
    cores         = max(1, min(4, parallel::detectCores(logical = FALSE) - 1, na.rm = TRUE)),
    batch_size    = NULL,                 # default: 3 x cores
    priority      = 15,                   # tools::psnice(): 15 = "below normal" on Windows, nice 15 elsewhere
    fdr_levels    = c(0.01, 0.05, 0.1),
    grid_filter   = NULL,                 # optional function(grid) -> logical, e.g. for a test
    max_batches   = Inf
  )
  user <- if (exists("sim2_config", envir = globalenv())) get("sim2_config", envir = globalenv()) else list()
  unknown <- setdiff(names(user), names(defaults))
  if (length(unknown)) stop("Unknown setting(s) in sim2_config: ", paste(unknown, collapse = ", "))
  cfg <- utils::modifyList(defaults, user, keep.null = TRUE)
  if (is.null(cfg$out_dir))     cfg$out_dir     <- file.path(cfg$pkg_dir, "benchmark", "results", "sim2")
  if (is.null(cfg$import_from)) cfg$import_from <- file.path(cfg$pkg_dir, "benchmark", "results")
  if (is.null(cfg$batch_size))  cfg$batch_size  <- 3 * cfg$cores
  cfg$cores <- max(1, as.integer(cfg$cores))
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  res_dir  <- dirname(cfg$out_dir)
  log_file <- file.path(cfg$out_dir, "progress.log")

  grid <- expand.grid(seed = seq_len(cfg$n_seeds), rho = c(0, 0.5, 0.9), pi_mcar = c(1, 0.5, 0.2),
                      rate = c(0.1, 0.2, 0.4), n_rep = c(3, 5), KEEP.OUT.ATTRS = FALSE)
  if (is.function(cfg$grid_filter)) grid <- grid[cfg$grid_filter(grid), ]
  grid$id <- sprintf("rep%d_rate%02d_mcar%03d_rho%02d_seed%02d", grid$n_rep, round(100 * grid$rate),
                     round(100 * grid$pi_mcar), round(10 * grid$rho), grid$seed)
  rds_of <- function(id) file.path(cfg$out_dir, paste0(id, ".rds"))

  # results of an earlier run saved in another folder (copied, never moved or deleted)
  if (dir.exists(cfg$import_from) && normalizePath(cfg$import_from) != normalizePath(cfg$out_dir)) {
    old <- file.path(cfg$import_from, paste0(grid$id, ".rds"))
    cp  <- file.exists(old) & !file.exists(rds_of(grid$id))
    if (any(cp)) {
      file.copy(old[cp], rds_of(grid$id[cp]))
      message(sum(cp), " finished datasets copied from ", cfg$import_from)
    }
  }

  # ---- one dataset (runs in a worker) -------------------------------------------------------
  run_one <- function(g, out_dir, msstats_seeds, fdr_levels, log_file) {
    f <- file.path(out_dir, paste0(g$id, ".rds"))
    if (file.exists(f)) return(f)
    t_start <- Sys.time()
    ok <- tryCatch({
      set.seed(g$seed)
      sim <- simulate_ptm_features(n_rep = g$n_rep, rho = g$rho, rate = c(g$rate / 2, g$rate), pi_mcar = g$pi_mcar, seed = g$seed)
      x0  <- filter_min_observed(sim_to_experiment(sim), min_n = 2)
      truth <- truth_table(x0)
      ctr <- make_contrasts(x0)
      have <- available_methods()
      arms <- list(); times <- list()
      run <- function(label, imputation, expr) {
        t0 <- Sys.time()
        r <- tryCatch(suppressWarnings(suppressMessages(expr)), error = function(e) { message(g$id, " ", label, ": ", conditionMessage(e)); NULL })
        if (!is.null(r) && nrow(r)) {
          r$method <- label; r$imputation <- imputation
          arms[[paste(label, imputation)]] <<- r
          times[[paste(label, imputation)]] <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        }
      }
      for (imp in c("none", "mle", "slsa")) {
        if (imp != "none" && !requireNamespace("imp4p", quietly = TRUE)) next
        x <- tryCatch(match_ptm_to_protein(impute(x0, method = imp), "id"), error = function(e) NULL)
        if (is.null(x)) next
        run("limma", imp, method_limma(x, ctr, pi0_method = "bh"))
        if ("msqrob2ptm" %in% have) run("msqrob2ptm", imp, method_msqrob2ptm(x, ctr))
        if (imp == "none") {
          run("limma_unpaired", imp, method_limma(x, ctr, analyses = "ptm_DPU", pi0_method = "bh", pair_by_sample = FALSE))
          run("limma_weighted", imp, method_limma(x, ctr, pi0_method = "bh", weights = TRUE))
          if (!is.null(arms[["limma none"]])) {
            nc <- arms[["limma none"]]; nc <- nc[nc$analysis == "ptm_DPA", ]; nc$analysis <- "ptm_DPU"
            nc$method <- "no_correction"; arms[["no_correction none"]] <- nc
          }
          if ("msstatsptm" %in% have && g$seed %in% msstats_seeds) run("msstatsptm", imp, method_msstatsptm(x, ctr))
        }
        rm(x)
      }
      if ("msstatsptm" %in% have && g$seed %in% msstats_seeds) {
        inp <- sim_to_msstats_input(sim)
        run("msstatsptm_native", "aft", method_msstatsptm_native(inp, ctr, MBimpute = TRUE))
        run("msstatsptm_native", "none", method_msstatsptm_native(inp, ctr, MBimpute = FALSE))
        rm(inp)
      }
      rm(sim, x0)
      res <- do.call(rbind, arms); rm(arms)
      res$arm <- paste(res$method, res$imputation, sep = "|")
      ev_res <- res; ev_res$method <- ev_res$arm                       # evaluate per arm
      split_arm <- function(d) { d$imputation <- sub("^.*\\|", "", d$method); d$method <- sub("\\|.*$", "", d$method); d }
      meta <- data.frame(id = g$id, seed = g$seed, n_rep = g$n_rep, rate = g$rate, pi_mcar = g$pi_mcar, rho = g$rho)
      out <- list(
        metrics  = cbind(meta, split_arm(evaluate_against_truth(ev_res, truth, fdr_levels))),
        by_class = cbind(meta, split_arm(evaluate_by_class(ev_res, truth, fdr_levels))),
        fdp      = cbind(meta, split_arm(fdp_curve(ev_res, truth))),
        runtime  = cbind(meta, data.frame(arm = names(times), seconds = unlist(times))))
      tmp <- paste0(f, ".part")                                      # write then rename: no half-written file
      saveRDS(out, tmp); file.rename(tmp, f)
      TRUE
    }, error = function(e) { message(g$id, ": ", conditionMessage(e)); FALSE })
    cat(sprintf("%s\t%s\t%s\t%.0f s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), g$id,
                if (ok) "done" else "FAILED", as.numeric(difftime(Sys.time(), t_start, units = "secs"))),
        file = log_file, append = TRUE)
    invisible(gc(verbose = FALSE))
    if (ok) f else NA_character_
  }

  # ---- summary tables, rewritten after each batch -------------------------------------------
  write_tables <- function() {
    files <- rds_of(grid$id); files <- files[file.exists(files)]
    if (!length(files)) return(0L)
    all <- lapply(files, readRDS)
    for (part in c("metrics", "by_class", "fdp", "runtime")) {
      name <- c(metrics = "sim2_metrics.tsv", by_class = "sim2_by_class.tsv", fdp = "sim2_fdp_curve.tsv",
                runtime = "sim2_runtime.tsv")[[part]]
      # datasets saved before 2026-11-05 name the method column "engine"
      tabs <- lapply(all, function(z) { d <- z[[part]]; if (!is.null(d)) names(d)[names(d) == "engine"] <- "method"; d })
      utils::write.table(do.call(rbind, tabs), file.path(res_dir, name),
                         sep = "\t", quote = FALSE, row.names = FALSE)
    }
    length(files)
  }

  init_worker <- function(pkg_dir, priority) {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    try(tools::psnice(value = priority), silent = TRUE)
    suppressMessages(devtools::load_all(pkg_dir, quiet = TRUE))
    if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)
    invisible(TRUE)
  }

  environment(run_one) <- globalenv(); environment(init_worker) <- globalenv()   # nothing large sent to workers

  # ---- run ----------------------------------------------------------------------------------
  todo <- grid[!file.exists(rds_of(grid$id)), ]
  n_total <- nrow(grid); n_done0 <- n_total - nrow(todo)
  message(sprintf("%d of %d datasets to run, %d worker(s), batches of %d. Progress: %s",
                  nrow(todo), n_total, cfg$cores, cfg$batch_size, log_file))
  if (nrow(todo)) {
    rows <- split(todo, seq_len(nrow(todo)))
    batches <- split(seq_along(rows), ceiling(seq_along(rows) / cfg$batch_size))
    t0 <- Sys.time(); n_run <- 0
    for (b in seq_along(batches)) {
      if (b > cfg$max_batches) { message("Stopped after ", cfg$max_batches, " batch(es); source the script again to resume."); break }
      these <- rows[batches[[b]]]
      if (cfg$cores > 1) {
        cl <- parallel::makeCluster(min(cfg$cores, length(these)))
        tryCatch({
          parallel::clusterCall(cl, init_worker, cfg$pkg_dir, cfg$priority)
          parallel::parLapplyLB(cl, these, run_one, out_dir = cfg$out_dir, msstats_seeds = cfg$msstats_seeds,
                                fdr_levels = cfg$fdr_levels, log_file = log_file)
        }, finally = parallel::stopCluster(cl))
      } else {
        if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)
        for (g in these) run_one(g, cfg$out_dir, cfg$msstats_seeds, cfg$fdr_levels, log_file)
      }
      n_run <- n_run + length(these)
      n_done <- write_tables()
      el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      left <- n_total - n_done
      msg <- sprintf("%s  batch %d/%d: %d/%d datasets done; %.0f min elapsed; about %.0f min left",
                     format(Sys.time(), "%H:%M"), b, length(batches), n_done, n_total, el,
                     el / n_run * left)
      message(msg); cat("#", msg, "\n", file = log_file, append = TRUE)
      gc(verbose = FALSE)
    }
  }
  n_done <- write_tables()
  message("Summary tables (", n_done, " of ", n_total, " datasets) written to ", res_dir)
  failed <- if (file.exists(log_file)) grep("FAILED", readLines(log_file), value = TRUE) else character(0)
  if (length(failed)) message(length(failed), " failed dataset(s) recorded in ", log_file,
                              " (they are retried at the next run)")
  invisible(NULL)
})

# Benchmark 2: MSstatsPTM ubiquitin spike-in (MassIVE MSV000088971, reanalysis RMSV000000669).
#
# Inputs (MassIVE.quant exports placed in benchmark/data/, see inst/benchmark/README.md):
#   summarised PTM intensities      -> table with columns RUN, Protein (site), LogIntensities, originalRUN, GROUP
#   summarised protein intensities  -> same columns, Protein = accession
#   MSstatsPTM adjusted model       -> columns Protein, Label, trueLog2FC, labeling (the ground truth)
# The three tables are located automatically from their headers.
#
# Truth: Kohler et al. 2023, Fig. 1B — heavy spiked sites: adjusted design log2FC; light sites: 0.
# Output: benchmark/results/spikein.tsv (+ spikein_by_contrast.tsv) and spikein_results.rds.

suppressPackageStartupMessages({ library(DiffPTM); library(data.table) })
data_dir <- "benchmark/data"; out_dir <- "benchmark/results"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
methods <- intersect(c("limma", "msstatsptm", "msqrob2ptm"), available_methods())

files <- list.files(data_dir, pattern = "display_quant_results-main\\.tsv$", recursive = TRUE, full.names = TRUE)
hdr <- lapply(files, function(f) names(fread(f, nrows = 0)))
is_summ <- vapply(hdr, function(h) all(c("originalRUN", "LogIntensities", "GROUP") %in% h), logical(1))
is_truth <- vapply(hdr, function(h) all(c("trueLog2FC", "protadj") %in% h) && !"issue" %in% h, logical(1))
summ <- lapply(files[is_summ], fread)
is_ptm <- vapply(summ, function(d) mean(grepl("_[A-Z]\\d+", d$Protein)) > 0.9, logical(1))
q <- summ[[which(is_ptm)[1]]]; p <- summ[[which(!is_ptm)[1]]]
tr <- fread(files[is_truth][1])
cat("PTM sites:", length(unique(q$Protein)), " proteins:", length(unique(p$Protein)), " truth rows:", nrow(tr), "\n")

to_matrix <- function(d) { w <- dcast(d, Protein ~ originalRUN, value.var = "LogIntensities", fun.aggregate = mean)
  m <- as.matrix(w[, -1]); rownames(m) <- w$Protein; m[is.nan(m)] <- NA; m }
P <- to_matrix(p); Q <- to_matrix(q)[, colnames(P)]
samples <- colnames(P); group <- sub("-.*$", "", samples)
ann <- make_annotation(samples, group)
prot <- new_assay(P, data.frame(feature_id = rownames(P), protein_ids = rownames(P), gene = NA_character_), log2_transform = FALSE)
site_prot <- sub("^([^_ ]+)_.*$", "\\1", rownames(Q))
ptm <- new_assay(Q, data.frame(feature_id = rownames(Q), protein_ids = site_prot, stripped_sequence = NA_character_,
                               modified_sequence = NA_character_, residue = sub("^[^_]+_([A-Z]).*$", "\\1", rownames(Q)),
                               position = sub("^[^_]+_[A-Z](\\d+).*$", "\\1", rownames(Q)), localization_prob = NA_real_,
                               ptm_type = "GlyGly", labeling = ifelse(grepl("heavy", rownames(Q)), "heavy", "light")),
                 log2_transform = FALSE)
x <- new_ptm_experiment(prot, ptm, annotation = list(protein = ann, ptm = ann)) |>
  filter_min_observed(min_n = 2) |> match_ptm_to_protein("id")      # already normalised by MSstats; no imputation
ctr <- make_contrasts(x, "pairwise")

# Truth (Kohler et al. 2023, Fig. 1): the 40 heavy spiked sites change by the design value of
# each contrast (0 for mix4 vs mix1); the light human sites are the background whose protein-
# adjusted change is 0 (human lysate is both the PTM background and the global proteome).
# The `trueLog2FC` column of the MassIVE table is the design value copied on every row, so
# it is only used for the heavy sites.
# Sign: our contrast "mixI_vs_mixJ" = mixI - mixJ ; MSstatsPTM label "mixJ vs mixI" = mixJ - mixI.
tr[, c("g_a", "g_b") := tstrsplit(Label, " vs ")]
tr[, true_adj := ifelse(labeling == "heavy", -trueLog2FC, 0)]
truth <- data.frame(analysis = "ptm_DPU", feature_id = tr$Protein, contrast = paste0(tr$g_b, "_vs_", tr$g_a),
                    positive = as.integer(tr$true_adj != 0), true_log2fc = tr$true_adj, stringsAsFactors = FALSE)
cat("truth: ", sum(truth$positive), "positive tests (heavy spiked sites),", sum(!truth$positive), "negatives\n")

res <- list()
for (e in methods) {
  t0 <- Sys.time()
  r <- tryCatch(suppressWarnings(run_methods(x, ctr, methods = e, pi0_method = "bh")), error = function(err) { message(e, ": ", conditionMessage(err)); NULL })
  if (is.null(r)) next
  r$runtime_s <- as.numeric(difftime(Sys.time(), t0, units = "secs")); res[[e]] <- r
  cat(sprintf("%-10s %6.1f s\n", e, r$runtime_s[1]))
}
# naive baseline: PTM abundance used as if it were usage
naive <- res$limma[res$limma$analysis == "ptm_DPA", ]; naive$analysis <- "ptm_DPU"; naive$method <- "no_correction"
# reference: the adjusted model published with the dataset (MSstatsPTM, Kohler et al. 2023)
pub <- data.frame(feature_id = tr$Protein, analysis = "ptm_DPU", contrast = paste0(tr$g_b, "_vs_", tr$g_a),
                  log2FC = -tr$log2FC, se = tr$SE, statistic = -tr$Tvalue, pvalue = tr$pvalue, adj_pvalue = tr$adj.pvalue,
                  method = "msstatsptm_published", runtime_s = NA_real_, stringsAsFactors = FALSE)
pub[!is.finite(pub$log2FC), c("log2FC", "pvalue", "adj_pvalue")] <- NA
res <- rbind(do.call(rbind, res), naive, pub)
saveRDS(list(x = x, res = res, truth = truth), file.path(out_dir, "spikein_results.rds"))

# Metrics. NOTE: after protein adjustment the light human sites are *not* a clean null on this
# dataset (every method, including the published model, calls 3-5 % of them), so
# `observed_fdr` is dominated by that background. Report sensitivity / log2FC accuracy on the
# heavy sites and the light-site call rate ("background_rate") side by side.
ev <- evaluate_against_truth(res, truth, fdr = c(0.01, 0.05))
ev_c <- evaluate_against_truth(res, truth, fdr = 0.05, by_contrast = TRUE)
m <- merge(res, truth, by = c("feature_id", "analysis", "contrast"))
bg <- aggregate(cbind(background_rate = adj_pvalue < 0.05) ~ method, m[m$positive == 0 & !is.na(m$adj_pvalue), ], mean)
ev <- merge(ev, bg, by = "method")
fwrite(ev, file.path(out_dir, "spikein.tsv"), sep = "\t"); fwrite(ev_c, file.path(out_dir, "spikein_by_contrast.tsv"), sep = "\t")
print(format(ev[order(ev$fdr, ev$method), ], digits = 3), row.names = FALSE)

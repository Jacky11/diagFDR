# Builds the small example exports in inst/extdata/formats/<software>/ used by the tests and
# the vignettes. One simulated experiment (2 conditions x 3 samples, 60 proteins, 45
# phosphosites) is written in the native layout of each software, with the column names of
# real exports (MaxQuant 2.x, Spectronaut 18/19, Proteome Discoverer 2.5/3.x, FragPipe 20-22,
# DIA-NN 1.8/1.9). Values are simulated: 10 proteins change (log2FC 1.5), 8 sites change in
# usage, 4 sites are compensated; ~10 % missing values.
# Usage (from the package folder): Rscript inst/scripts/make_format_fixtures.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(42)
out_root <- "inst/extdata/formats"
n_prot <- 60; n_site <- 45; cond <- rep(c("control", "treated"), each = 3); rep_id <- rep(1:3, 2)
samples <- paste0(cond, "_", rep_id)
acc <- sprintf("P%05d", sample(10000:99999, n_prot)); gene <- paste0("GENE", seq_len(n_prot))
aa <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
prot_seq <- vapply(seq_len(n_prot), function(i) paste(sample(aa, 180, TRUE), collapse = ""), "")
# proteins
base_p <- rnorm(n_prot, 25, 1.5)
P <- matrix(base_p, n_prot, 6) + matrix(rnorm(n_prot * 6, 0, 0.25), n_prot)
P[1:10, cond == "treated"] <- P[1:10, cond == "treated"] + 1.5
dimnames(P) <- list(acc, samples)
# sites: site k belongs to protein k (k <= 45); tryptic-like peptide around a random S/T/Y
site_prot <- seq_len(n_site)
pep <- character(n_site); pos <- integer(n_site); res <- character(n_site); pep_start <- integer(n_site); pos_in_pep <- integer(n_site)
for (k in seq_len(n_site)) {
  s <- prot_seq[site_prot[k]]
  st <- sample(20:150, 1); len <- sample(8:14, 1)
  p <- substr(s, st, st + len - 1); r <- sample(c("S", "T", "Y"), 1); ip <- sample(2:(len - 1), 1)
  substr(p, ip, ip) <- r; substr(s, st + ip - 1, st + ip - 1) <- r
  substr(p, len, len) <- "K"; substr(s, st + len - 1, st + len - 1) <- "K"; substr(s, st - 1, st - 1) <- "R"
  prot_seq[site_prot[k]] <- s; pep[k] <- p; pep_start[k] <- st; pos_in_pep[k] <- ip; pos[k] <- st + ip - 1; res[k] <- r
}
Q <- P[site_prot, ] - 3 + matrix(rnorm(n_site * 6, 0, 0.3), n_site)
Q[1:4, cond == "treated"] <- Q[1:4, cond == "treated"] - 1.5      # compensated: protein up, site flat -> usage down
Q[11:18, cond == "treated"] <- Q[11:18, cond == "treated"] + 1.5  # usage up
dimnames(Q) <- list(paste0(acc[site_prot], "_", res, pos), samples)
miss <- function(M, rate = 0.1) { M[sample(length(M), round(rate * length(M)))] <- NA; M }
P <- miss(P); Q <- miss(Q, 0.15)
Pi <- round(2^P); Qi <- round(2^Q)           # raw intensities
flank <- vapply(seq_len(n_site), function(k) { s <- prot_seq[site_prot[k]]; w <- substr(paste0(strrep("_", 7), s, strrep("_", 7)), pos[k], pos[k] + 14); w }, "")
fasta <- unlist(lapply(seq_len(n_prot), function(i) c(sprintf(">sp|%s|%s_HUMAN Protein %d OS=Homo sapiens GN=%s", acc[i], gene[i], i, gene[i]), prot_seq[i])))
dir.create(out_root, showWarnings = FALSE, recursive = TRUE); writeLines(fasta, file.path(out_root, "example.fasta"))
na0 <- function(M) { M[is.na(M)] <- 0; M }
wtsv <- function(d, f) utils::write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# ---- MaxQuant: combined search, proteome runs + enriched runs (Phospho_ prefix) ------------
d <- file.path(out_root, "MaxQuant"); dir.create(d, showWarnings = FALSE)
pg <- data.frame(check.names = FALSE, `Protein IDs` = acc, `Majority protein IDs` = acc, `Protein names` = paste("Protein", seq_len(n_prot)),
                 `Gene names` = gene, `Fasta headers` = fasta[seq(1, length(fasta), 2)], `Number of proteins` = 1, Peptides = sample(3:30, n_prot, TRUE),
                 `Razor + unique peptides` = sample(3:30, n_prot, TRUE), `Unique peptides` = sample(1:20, n_prot, TRUE), `Sequence coverage [%]` = round(runif(n_prot, 10, 80), 1),
                 `Mol. weight [kDa]` = round(runif(n_prot, 20, 150), 2), `Q-value` = 0, Score = round(runif(n_prot, 20, 300), 2), Intensity = rowSums(na0(Pi)))
for (s in samples) pg[[paste("Intensity", s)]] <- na0(Pi)[, s]
for (s in samples) pg[[paste("Intensity", paste0("Phospho_", s))]] <- round(na0(Pi)[, s] * runif(n_prot, 0, 0.05))
for (s in samples) pg[[paste("LFQ intensity", s)]] <- na0(Pi)[, s]
pg$`MS/MS count` <- sample(5:100, n_prot, TRUE); pg$id <- seq_len(n_prot) - 1; pg$`Only identified by site` <- ""; pg$Reverse <- ""; pg$`Potential contaminant` <- ""
pg <- rbind(pg, within(pg[1:2, ], { `Protein IDs` <- c("REV__P99999", "CON__P88888"); `Majority protein IDs` <- `Protein IDs`; Reverse <- c("+", ""); `Potential contaminant` <- c("", "+") }))
wtsv(pg, file.path(d, "proteinGroups.txt"))
st <- data.frame(check.names = FALSE, Proteins = acc[site_prot], `Positions within proteins` = pos, `Leading proteins` = acc[site_prot], Protein = acc[site_prot],
                 `Protein names` = paste("Protein", site_prot), `Gene names` = gene[site_prot], `Fasta headers` = fasta[2 * site_prot - 1],
                 `Localization prob` = round(c(runif(40, 0.8, 1), runif(5, 0.3, 0.7)), 4), PEP = signif(runif(n_site, 1e-8, 1e-3), 3), Score = round(runif(n_site, 50, 200), 2),
                 `Delta score` = round(runif(n_site, 5, 100), 2), `Score for localization` = round(runif(n_site, 50, 200), 2))
for (s in samples) st[[paste("Localization prob", paste0("Phospho_", s))]] <- st$`Localization prob`
st$Position <- pos; st$`Amino acid` <- res; st$`Sequence window` <- flank; st$`Modification window` <- "X"; st$`Peptide window coverage` <- 1
st$`Phospho (STY) Probabilities` <- pep; st$`Position in peptide` <- pos_in_pep; st$Charge <- 2; st$`Mass error [ppm]` <- round(rnorm(n_site, 0, 0.5), 3)
st$Intensity <- rowSums(na0(Qi)); st$Intensity___1 <- st$Intensity; st$Intensity___2 <- 0; st$Intensity___3 <- 0; st$`Ratio mod/base` <- NA
for (s in samples) st[[paste("Intensity", s)]] <- 0
for (s in samples) st[[paste("Intensity", paste0("Phospho_", s))]] <- na0(Qi)[, s]
for (s in samples) for (m in 1:3) st[[paste0("Intensity ", s, "___", m)]] <- 0
for (s in samples) for (m in 1:3) st[[paste0("Intensity Phospho_", s, "___", m)]] <- if (m == 1) na0(Qi)[, s] else if (m == 2) round(na0(Qi)[, s] * ifelse(seq_len(n_site) %% 5 == 0, 0.3, 0)) else 0
st$Reverse <- ""; st$`Potential contaminant` <- ""; st$id <- seq_len(n_site) - 1; st$`Protein group IDs` <- site_prot - 1; st$Positions <- pos; st$`Peptide IDs` <- seq_len(n_site)
wtsv(st, file.path(d, "Phospho_STY_Sites.txt"))

# ---- Spectronaut: semicolon-separated pivot reports with decimal commas ------------------------
d <- file.path(out_root, "Spectronaut"); dir.create(d, showWarnings = FALSE)
fmt <- function(v) ifelse(is.na(v), "Filtered", sub(".", ",", format(v, scientific = FALSE, trim = TRUE), fixed = TRUE))
run_p <- sprintf("[%d] Proteome-%02d_%s_R%d", seq_along(samples), seq_along(samples), cond, rep_id)
run_s <- sprintf("[%d] Phospho-%02d_%s_R%d", seq_along(samples), seq_along(samples), cond, rep_id)
grp <- acc; grp[3] <- paste0(acc[3], ";", "Q99999")          # one protein group with two accessions (quoted, on two lines in the file)
pr <- data.frame(check.names = FALSE, PG.ProteinGroups = grp, PG.Genes = gene, PG.ProteinDescriptions = paste("Protein", seq_len(n_prot)),
                 PG.FastaHeaders = fasta[seq(1, length(fasta), 2)], PG.Qvalue = fmt(signif(runif(n_prot, 1e-6, 1e-3), 3)), PG.MolecularWeight = fmt(round(runif(n_prot, 20, 150), 2)))
for (i in seq_along(samples)) pr[[paste0(run_p[i], ".PG.NrOfStrippedSequencesIdentified")]] <- sample(2:30, n_prot, TRUE)
for (i in seq_along(samples)) pr[[paste0(run_p[i], ".PG.Quantity")]] <- fmt(round(2^P[, i] / 10, 2))
for (i in seq_along(samples)) pr[[paste0(run_p[i], ".PG.IBAQ")]] <- fmt(round(2^P[, i] / 1000, 3))
lines <- c(paste(names(pr), collapse = ";"), apply(pr, 1, function(r) { r <- trimws(r); paste(ifelse(grepl(";", r), paste0('"', sub(";", "\n", r), '"'), r), collapse = ";") }))
writeLines(lines, file.path(d, "PivotReport_Proteomes.csv"))
sp_rows <- lapply(seq_len(n_site), function(k) {
  accs <- if (k == 3) c(acc[site_prot[k]], "Q99999") else acc[site_prot[k]]       # site of a protein group: one row per accession
  lapply(accs, function(a) data.frame(check.names = FALSE, PG.ProteinGroups = grp[site_prot[k]], PG.Genes = gene[site_prot[k]], PTM.ProteinId = a,
    PTM.CollapseKey = sprintf("%s_%s%d_M1", a, res[k], pos[k]), PTM.Multiplicity = 1, PTM.ModificationTitle = "Phospho (STY)", PTM.SiteAA = res[k],
    PTM.SiteLocation = pos[k], PTM.FlankingRegion = flank[k], stringsAsFactors = FALSE))
})
sp <- do.call(rbind, unlist(sp_rows, recursive = FALSE)); kk <- rep(seq_len(n_site), vapply(sp_rows, length, 1L))
for (i in seq_along(samples)) sp[[paste0(run_s[i], ".PTM.Group")]] <- ifelse(is.na(Q[kk, i]), "Filtered", paste0("_", pep[kk], "_.2"))
for (i in seq_along(samples)) sp[[paste0(run_s[i], ".PTM.SiteProbability")]] <- fmt(ifelse(is.na(Q[kk, i]), NA, st$`Localization prob`[kk]))
for (i in seq_along(samples)) sp[[paste0(run_s[i], ".PTM.Quantity")]] <- fmt(round(2^Q[kk, i] / 10, 2))
lines <- c(paste(names(sp), collapse = ";"), apply(sp, 1, function(r) { r <- trimws(r); paste(ifelse(grepl(";", r), paste0('"', sub(";", "\n", r), '"'), r), collapse = ";") }))
writeLines(lines, file.path(d, "PivotReport_Phospho.csv"))

# ---- Proteome Discoverer: Proteins and PeptideGroups exports ---------------------------------------
d <- file.path(out_root, "ProteomeDiscoverer"); dir.create(d, showWarnings = FALSE)
pd_abund <- function(M, prefix) { out <- list(); for (i in seq_along(samples)) out[[sprintf("%s: F%d: Sample, %s", prefix, i, cond[i])]] <- round(M[, i], 1); as.data.frame(out, check.names = FALSE) }
pdp <- data.frame(check.names = FALSE, `Protein FDR Confidence: Combined` = "High", Master = "IsMasterProtein", Accession = acc, Description = sprintf("Protein %d OS=Homo sapiens GN=%s", seq_len(n_prot), gene),
                  `Exp. q-value: Combined` = 0, `Sum PEP Score` = round(runif(n_prot, 10, 300), 3), Contaminant = "False", `Coverage [%]` = round(runif(n_prot, 10, 80), 0),
                  `# Peptides` = sample(3:30, n_prot, TRUE), `# PSMs` = sample(5:100, n_prot, TRUE), `# Unique Peptides` = sample(1:20, n_prot, TRUE), `# Protein Groups` = 1,
                  `# AAs` = 180, `MW [kDa]` = round(runif(n_prot, 20, 150), 1), `calc. pI` = round(runif(n_prot, 4, 10), 2), `Gene Symbol` = gene)
pdp <- cbind(pdp, pd_abund(2^P, "Abundance"), pd_abund(2^P * 1.02, "Abundances (Normalized)"))
pdp$`Abundances (Grouped): control` <- round(rowMeans(2^P[, 1:3], na.rm = TRUE), 1); pdp$`Abundances (Grouped): treated` <- round(rowMeans(2^P[, 4:6], na.rm = TRUE), 1)
pdp$`Abundance Ratio: (treated) / (control)` <- round(pdp$`Abundances (Grouped): treated` / pdp$`Abundances (Grouped): control`, 3)
wtsv(pdp, file.path(d, "Proteins.txt"))
pdq <- data.frame(check.names = FALSE, Checked = "False", Confidence = "High", `Annotated Sequence` = sprintf("[R].%s.[A]", pep), Sequence = pep,
                  Modifications = sprintf("1xPhospho [%s%d(%s)]", res, pos_in_pep, round(100 * st$`Localization prob`, 1)),
                  `Modifications in Master Proteins` = sprintf("%s 1xPhospho [%s%d(%s)]", acc[site_prot], res, pos, round(100 * st$`Localization prob`, 1)),
                  `Qvality PEP` = signif(runif(n_site, 1e-6, 1e-2), 3), `Qvality q-value` = 0, `# Protein Groups` = 1, `# Proteins` = 1, `# PSMs` = sample(1:10, n_site, TRUE),
                  `Master Protein Accessions` = acc[site_prot], `Master Protein Descriptions` = paste("Protein", site_prot),
                  `Positions in Master Proteins` = sprintf("%s [%d-%d]", acc[site_prot], pep_start, pep_start + nchar(pep) - 1), `# Missed Cleavages` = 0, Contaminant = "False")
pdq <- cbind(pdq, pd_abund(2^Q, "Abundance"), pd_abund(2^Q * 1.02, "Abundances (Normalized)"))
wtsv(pdq, file.path(d, "PeptideGroups.txt"))

# ---- FragPipe label-free: combined_protein.tsv + combined_modified_peptide.tsv ------------------
d <- file.path(out_root, "FragPipe_LFQ"); dir.create(d, showWarnings = FALSE)
fp <- data.frame(check.names = FALSE, Protein = sprintf("sp|%s|%s_HUMAN", acc, gene), `Protein ID` = acc, `Entry Name` = paste0(gene, "_HUMAN"), Gene = gene, `Protein Length` = 180,
                 Organism = "Homo sapiens OX=9606", `Protein Existence` = "1:Experimental evidence at protein level", Description = paste("Protein", seq_len(n_prot)),
                 `Protein Probability` = 1, `Top Peptide Probability` = 1, `Combined Total Peptides` = sample(3:30, n_prot, TRUE), `Combined Spectral Count` = sample(5:100, n_prot, TRUE),
                 `Combined Unique Spectral Count` = sample(5:100, n_prot, TRUE), `Combined Total Spectral Count` = sample(5:100, n_prot, TRUE))
for (s in samples) { fp[[paste(s, "Spectral Count")]] <- sample(0:30, n_prot, TRUE); fp[[paste(s, "Unique Spectral Count")]] <- sample(0:30, n_prot, TRUE); fp[[paste(s, "Total Spectral Count")]] <- sample(0:30, n_prot, TRUE) }
for (s in samples) fp[[paste(s, "Intensity")]] <- round(na0(Pi)[, s] * 1.1)
for (s in samples) fp[[paste(s, "MaxLFQ Intensity")]] <- na0(Pi)[, s]
fp$`Indistinguishable Proteins` <- ""
wtsv(fp, file.path(d, "combined_protein.tsv"))
modseq <- vapply(seq_len(n_site), function(k) { p <- pep[k]; paste0(substr(p, 1, pos_in_pep[k]), "[", 79.9663 + ifelse(res[k] == "S", 87.0320, ifelse(res[k] == "T", 101.0477, 163.0633)), "]", substr(p, pos_in_pep[k] + 1, nchar(p))) }, "")
fq <- data.frame(check.names = FALSE, `Peptide Sequence` = pep, `Modified Sequence` = modseq, `Prev AA` = "R", `Next AA` = "A", Start = pep_start, End = pep_start + nchar(pep) - 1,
                 `Peptide Length` = nchar(pep), Charges = "2, 3", Probability = round(st$`Localization prob`, 4), `Spectral Count` = sample(1:20, n_site, TRUE), Intensity = rowSums(na0(Qi)),
                 `Assigned Modifications` = sprintf("%d%s(79.9663)", pos_in_pep, res), `Observed Modifications` = "", Protein = sprintf("sp|%s|%s_HUMAN", acc[site_prot], gene[site_prot]),
                 `Protein ID` = acc[site_prot], `Entry Name` = paste0(gene[site_prot], "_HUMAN"), Gene = gene[site_prot], `Protein Description` = paste("Protein", site_prot), `Mapped Genes` = "", `Mapped Proteins` = "")
for (s in samples) fq[[paste(s, "Spectral Count")]] <- sample(0:10, n_site, TRUE)
for (s in samples) fq[[paste(s, "Intensity")]] <- round(na0(Qi)[, s] * 1.1)
for (s in samples) fq[[paste(s, "MaxLFQ Intensity")]] <- na0(Qi)[, s]
wtsv(fq, file.path(d, "combined_modified_peptide.tsv"))

# ---- FragPipe TMT (TMT-Integrator): abundance_protein_MD.tsv + abundance_single-site_MD.tsv (log2 ratios to reference) ----
d <- file.path(out_root, "FragPipe_TMT"); dir.create(d, showWarnings = FALSE)
refp <- rowMeans(P, na.rm = TRUE); refq <- rowMeans(Q, na.rm = TRUE)
tp <- data.frame(check.names = FALSE, Index = acc, NumberPSM = sample(5:100, n_prot, TRUE), ProteinID = acc, MaxPepProb = 1, ReferenceIntensity = round(refp, 4))
for (i in seq_along(samples)) tp[[samples[i]]] <- round(P[, i], 4)
wtsv(tp, file.path(d, "abundance_protein_MD.tsv"))
tq <- data.frame(check.names = FALSE, Index = sprintf("%s_%s%d", acc[site_prot], res, pos), Gene = gene[site_prot], ProteinID = acc[site_prot], Peptide = pep,
                 SequenceWindow = flank, Start = pep_start, End = pep_start + nchar(pep) - 1, MaxPepProb = round(st$`Localization prob`, 4), ReferenceIntensity = round(refq, 4))
for (i in seq_along(samples)) tq[[samples[i]]] <- round(Q[, i], 4)
wtsv(tq, file.path(d, "abundance_single-site_MD.tsv"))

# ---- DIA-NN: report.pg_matrix.tsv + report.pr_matrix.tsv (columns named by raw file paths) ------
d <- file.path(out_root, "DIA-NN"); dir.create(d, showWarnings = FALSE)
raw <- sprintf("D:\\\\data\\\\run_%s.raw", samples)
dg <- data.frame(check.names = FALSE, Protein.Group = acc, Protein.Ids = acc, Protein.Names = paste0(gene, "_HUMAN"), Genes = gene, First.Protein.Description = paste("Protein", seq_len(n_prot)))
for (i in seq_along(samples)) dg[[raw[i]]] <- round(2^P[, i], 1)
wtsv(dg, file.path(d, "report.pg_matrix.tsv"))
dmod <- vapply(seq_len(n_site), function(k) paste0(substr(pep[k], 1, pos_in_pep[k]), "(UniMod:21)", substr(pep[k], pos_in_pep[k] + 1, nchar(pep[k]))), "")
dq <- data.frame(check.names = FALSE, Protein.Group = acc[site_prot], Protein.Ids = acc[site_prot], Protein.Names = paste0(gene[site_prot], "_HUMAN"), Genes = gene[site_prot],
                 First.Protein.Description = paste("Protein", site_prot), Proteotypic = 1, Stripped.Sequence = pep, Modified.Sequence = dmod, Precursor.Charge = 2, Precursor.Id = paste0(dmod, "2"))
for (i in seq_along(samples)) dq[[raw[i]]] <- round(2^Q[, i], 1)
# a second charge state for the first ten peptidoforms (DIA-NN reports precursors)
dq2 <- dq[1:10, ]; dq2$Precursor.Charge <- 3; dq2$Precursor.Id <- paste0(dmod[1:10], "3"); for (i in seq_along(samples)) dq2[[raw[i]]] <- round(2^Q[1:10, i] * 0.4, 1)
wtsv(rbind(dq, dq2), file.path(d, "report.pr_matrix.tsv"))

# ---- annotation shared by all formats ----------------------------------------------------------------
utils::write.table(data.frame(sample = samples, condition = cond, subject = rep_id), file.path(out_root, "annotation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(protein_changed = acc[1:10], usage_up = rownames(Q)[11:18], compensated = rownames(Q)[1:4]), file.path(out_root, "truth.rds"))
cat("written", out_root, "\n"); print(system(paste("du -sh", out_root), intern = TRUE))

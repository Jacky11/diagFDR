# DiffPTM — help

## What the app does

A modified peptide's intensity changes when the modification changes **and** when its protein
changes. DiffPTM separates the two:

* **Differential PTM abundance (DPA)** — does the site (or modified peptide) change between
  conditions, whatever the reason?
* **Differential PTM usage (DPU)** — does the site change *relative to its parent protein*
  (is the stoichiometry of the modification altered)? This is the test that needs the
  total-proteome table.

Three statistical methods answer the DPU question and can be run side by side on the same
processed data:

| Method | How | Reference |
|---|---|---|
| Single model (limma) | site and protein intensities fitted together; the condition × component interaction is the usage change. Paired (uses the within-sample site–protein correlation) when PTM and proteome samples share the same names | Giai Gianetto, Methods Mol Biol 2023; Wolski et al., submitted; Ritchie et al. 2015 |
| MSstatsPTM | site and protein modelled separately, then the protein fold change is subtracted and the uncertainties added | Kohler et al., MCP 2023 |
| msqrob2PTM | protein intensity subtracted from the site in each sample, then the differences are modelled (robust regression with variance moderation) | Demeulemeester et al., MCP 2024 |

Estimates of the usage change are the same for the three; their standard errors differ (pairing,
variance moderation), hence different numbers of significant sites. With few replicates and
paired samples the single model is the most sensitive; MSstatsPTM is the choice when the PTM
and proteome samples are not the same.

## The six steps

1. **Import** — select the software and upload the proteome and PTM tables (see formats
   below). Tick off "Intensities are raw" for tables already on the log2 scale (FragPipe
   TMT-Integrator); the app warns when the scale looks wrong.
2. **Design** — DiffPTM proposes the design from the data: runs without any value are not
   used, enriched and total runs of a combined search are told apart by their number of sites,
   PTM and proteome runs of the same biological sample get the same *sample* name (the parts of
   the names that only mark the fraction, such as `Phospho_`, `P_`, `Proteome-`, are removed;
   these samples are paired in the usage test), and conditions are read from the sample names
   without the replicate number. When the conditions combine two factors (e.g. mock / SARS_COV2 x
   6h / 24h), choose whether to compare the conditions within each time point (default), the other
   way round, or all groups as one factor. Check the summary, edit cells if needed (double-click;
   an empty condition means the run is not used) and click *Apply design*. *Advanced options*:
   bulk renaming, other parts of the names as condition, second factor, subject, TMT mixture,
   or an annotation file.
3. **Processing** — remove reverse hits / contaminants, merge site rows repeated for
   several protein accessions (Spectronaut), filter on localization probability and on
   observed values, normalise (global median centering by default), match sites to proteins
   (by identifier, or by locating the peptide in the FASTA, which also gives site positions for
   peptide tables), optionally collapse peptidoforms to sites. Each option is explained on the
   page; above the normalisation options, DiffPTM says what the data suggest
   (`normalisation_diagnostics()`): a shift of the PTM assay between conditions confined to a
   subset of sites calls for per-condition centering; no shift, or a shift of every site
   (loading, enrichment efficiency, or a global biological change that the data cannot tell
   apart), for global centering. Imputation is optional and off by default: the tests use the observed values. The
   quality-control plots show the number of features observed per sample (failed runs stand
   out), the intensity distributions, a PCA and the variance of its components, the correlation
   matrix and the Jaccard index of detection between samples with the clustering of the samples
   on each, and the overlap of the features detected in the samples (UpSet plot for all samples,
   Venn diagram for the samples of a condition; features or modified proteins).
4. **Statistics** — choose the methods, the analyses, the contrasts (each condition against a
   reference, every pair, or each condition against the others; within each level of the
   second factor if there is one), the FDR and an optional fold-change threshold, and click
   *Run*. MSstatsPTM is much slower than the two other forms (11 min vs 12 s and 1 min for
   21,770 sites and four contrasts).
5. **Results** — volcano plot with gene/site labels, p-value histogram, PTM-vs-protein fold
   changes (points off the diagonal change in usage), counts per contrast and method,
   agreement between methods (UpSet plot), the table of significant features, and the PTM
   features detected in samples where their protein is not (detection discordance).
6. **Download** — Excel workbook (all tables, a README sheet with the Materials & Methods
   paragraph, which can be edited on the page first, and a sheet explaining every column),
   figures as editable PowerPoint and PDF (quality control first), and the session as an `.rds`
   file that `readRDS()` reloads in R.

## Input formats

| Software | Proteome table | PTM table | Notes |
|---|---|---|---|
| MaxQuant | `proteinGroups.txt` (LFQ intensities when present) | `Phospho (STY)Sites.txt`, `GlyGly (K)Sites.txt`, ... | combined searches list the proteome runs in the site table too: rename the `Phospho_` runs in step 2; multiplicity columns `___1/2/3` are expanded, kept or summed |
| Spectronaut | protein pivot report (`PG.Quantity`) | PTM site pivot report (`PTM.Quantity`, `PTM.SiteProbability`, `PTM.SiteAA`, `PTM.SiteLocation`) | semicolon separators and decimal commas are detected; `Filtered` = missing; one row per accession of a protein group is merged in step 3 |
| Proteome Discoverer | Proteins export | PeptideGroups export | tab-separated text exports; per-file abundances (`Abundance: F1: Sample`), normalised ones preferred; sites read from "Modifications in Master Proteins" |
| FragPipe | `combined_protein.tsv` (label-free) or `abundance_protein_*.tsv` (TMT) | `combined_modified_peptide.tsv` (label-free) or `abundance_single-site_*.tsv` (TMT) | TMT-Integrator tables are log2 ratios: untick "Intensities are raw" |
| DIA-NN | `report.pg_matrix.tsv` | `report.pr_matrix.tsv` (precursors; collapse to sites in step 3 with the FASTA) or a site matrix | columns are named after the raw files; the `.raw` / `.d` extension is removed |
| Generic | wide table: `protein_id`, `gene`, `intensity_<sample>` | wide table: `feature_id`, `protein_ids`, `gene`, `stripped_sequence`, `modified_sequence`, `residue`, `position`, `localization_prob`, `ptm_type`, `intensity_<sample>` | any other software, after reshaping |

"Advanced: column mapping" accepts a regular expression with one capture group naming the
sample, for exports whose column names differ from the defaults. The package vignettes
(`browseVignettes("DiffPTM")`) walk through one dataset per software.

## How to cite

If you use DiffPTM, please cite:

* the package: Giai Gianetto Q, Salhi Z. *DiffPTM: differential abundance and usage analysis of
  post-translational modifications.* R package, https://github.com/Jacky11/DiffPTM (in R:
  `citation("DiffPTM")`);
* the workflow it extends (quality control, detection profiles of PTMs relative to their protein,
  and the single model testing whether a PTM evolves differently from its protein across several
  biological conditions, illustrated on SARS-CoV-2 infected cells): Giai Gianetto Q. Statistical
  Analysis of Post-Translational Modifications Quantified by Label-Free Proteomics Across Multiple
  Biological Conditions with R: Illustration from SARS-CoV-2 Infected Cells. In: Burger T (ed.),
  *Statistical Analysis of Proteomic Data*. Methods Mol Biol. 2023;2426:267–302.
  [doi:10.1007/978-1-0716-1967-4_12](https://doi.org/10.1007/978-1-0716-1967-4_12);
* the model forms you used (table above), and cp4p (Giai Gianetto et al., Proteomics 2016) for the
  FDR control.

The Materials & Methods text of the Download step includes these references.

## Limitations

* MSstatsPTM and msqrob2PTM are optional dependencies; when they are not installed the
  corresponding options are disabled. MSstatsPTM takes minutes for thousands of sites.
* Sequence-based matching and site positions of peptide tables require the FASTA used for the
  searches.
* Multisite features are grouped by the set of sites only (no probability weighting).

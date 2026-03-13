<h1 align="center">diagFDR: Verifiable FDR Diagnostics for Proteomics</h1>

<div align="center" style="margin-top: 10px;">
  <a href="https://www.gnu.org/licenses/gpl-3.0"><img src="https://img.shields.io/badge/License-GPLv3-blue.svg" alt="GPL v3 License"></a>
</div></a>
  
## 1. Description

This package provides R functions to compute verifiable false discovery rate (FDR) diagnostic checks for workflows based on target-decoy competition and related confidence measures.
Implements calibration, stability and tail diagnostics, including tail support, threshold elasticity, posterior error probability (PEP) reliability, and equal-chance checks.

## 2. Installation

It is recommended to install the latest version of R. The installation of the `diagFDR` package can be performed by:

```
if(!requireNamespace("devtools", quietly = TRUE)){
   install.packages("devtools") 
}
devtools::install_github("Jacky11/diagFDR")
```
You can now load it and run the app with this commands:
```
library(diagFDR)
```

## 3. Verifiable FDR diagnostics HTML report

An standardized HTML report can be obtained from any identification software using the dfdr_render_report function to check the scope, the calibration, and the stability of a FDR-based identification search. Here we provide the way to get it from most common software.

### 3.1 DIA-NN

To enable all diagnostics, export:

- decoys: `--report-decoys`
- a permissive export ceiling: `--qvalue 0.5` (or higher)

The q-value ceiling matters because some diagnostics operate in low-confidence regions
(e.g. equal-chance plausibility checks, or local-window support around cutoffs).

First, you have to upload the .parquet file in your R session:
```
rep <- read_diann_parquet("path/to/report.parquet")
```

# (A) Global precursor list using Global.Q.Value
# Recommended for experiment-wide (pooled) lists.
x_global_gq <- diann_global_precursor(
  rep,
  q_col = "Global.Q.Value",
  q_max_export = 0.5,
  unit = "precursor",
  scope = "global",
  q_source = "Global.Q.Value"
)

# (B) Run×precursor universe using run-wise Q.Value
# Recommended for per-run decisions / QC.
x_runx <- diann_runxprecursor(
  rep,
  q_col = "Q.Value",
  q_max_export = 0.5,
  id_mode = "runxid",
  unit = "runxprecursor",
  scope = "runwise",
  q_source = "Q.Value"
)

# (C) Scope misuse comparator: min run-wise q over runs per precursor (anti-pattern)
# Useful for demonstrating/diagnosing scope mismatch.
x_minrun <- diann_global_minrunq(
  rep,
  q_col = "Q.Value",
  q_max_export = 0.5,
  unit = "precursor",
  scope = "aggregated",
  q_source = "min_run(Q.Value)"
)

diag <- dfdr_run_all(
  xs = list(global = x_global_gq, runx = x_runx, minrun = x_minrun),
  alpha_main = 0.01,
  compute_pseudo_pvalues = TRUE  # <-- This adds p-value diagnostics
)

# Compare accepted lists across scopes (Jaccard overlap across alpha)
scope_tbl <- dfdr_scope_disagreement(
  x1 = x_global_gq,
  x2 = x_minrun,
  alphas = c(1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2),
  label1 = "Global.Q.Value",
  label2 = "min_run(Q.Value)"
)

# Write outputs to disk (tables + plots; optionally PPTX)
dfdr_write_report(diag, out_dir = "diagFDR_diann_out", formats = c("csv", "png", "manifest", "readme", "summary"))

# Render a single HTML report (requires rmarkdown in Suggests)
dfdr_render_report(diag, out_dir = "diagFDR_diann_out")

---

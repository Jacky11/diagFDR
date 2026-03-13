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

## 3. Verifiable FDR diagnostics HTML reports

A standardized HTML report can be obtained from any identification software using the dfdr_render_report() function to check the scope, the calibration, and the stability of a FDR-based identification search. Here we provide the way to get it from most common software.

### 3.1 DIA-NN

To enable all diagnostics at the precursor level, perform the search in DIA-NN using:

- decoys export: `--report-decoys`
- a permissive export ceiling: `--qvalue 0.5` (or higher)

The q-value ceiling matters because some diagnostics operate in low-confidence regions
(e.g. equal-chance plausibility checks, or local-window support around cutoffs).

First, you have to upload the  `.parquet ` file in your R session:
```
rep <- read_diann_parquet("path/to/report.parquet")
```

Next, different kinds of scopes and q-values are available. 
A global precursor universe using `Global.Q.Value` recommended for experiment-wide (pooled) lists can be extracted:
```
x_global_gq <- diann_global_precursor(
  rep,
  q_col = "Global.Q.Value",
  q_max_export = 0.5,
  unit = "precursor",
  scope = "global",
  q_source = "Global.Q.Value"
)
```
A run×precursor universe using run-wise Q.Value recommended for per-run decisions can also be extracted:
```
x_runx <- diann_runxprecursor(
  rep,
  q_col = "Q.Value",
  q_max_export = 0.5,
  id_mode = "runxid",
  unit = "runxprecursor",
  scope = "runwise",
  q_source = "Q.Value"
)
```
All diagnostics measures can be obtained with:
```
diag <- dfdr_run_all(
  xs = list(global = x_global_gq, runx = x_runx, minrun = x_minrun),
  alpha_main = 0.01,
  compute_pseudo_pvalues = TRUE  # <-- This adds p-value diagnostics
)
```
The lists can be compared across scopes (Jaccard overlap across alpha) using:
```
scope_tbl <- dfdr_scope_disagreement(
  x1 = x_global_gq,
  x2 = x_runx,
  alphas = c(1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2),
  label1 = "Global.Q.Value",
  label2 = "Q.Value"
)
```

Finally, a single HTML report can be obtained using:
```
dfdr_render_report(diag, out_dir = "path/to/diagFDR_diann_out")
```

### 3.2 MaxQuant

To enable all diagnostics at the precursor level, the `msms.txt` file resulting from the search has to be used.

First, you have to upload the  `msms.txt ` file in your R session:
```
mq_path <- "path/to/msms.txt"

# Read msms.txt and reconstruct TDC q-values using MaxQuant Score and Reverse indicator.
# - Reverse == "+" is treated as a decoy indicator.
# - Score is assumed "higher is better".
# - q-values are computed using FDR(i) = (D(i)+add_decoy)/T(i) and q(i)=min_{j>=i} FDR(j).
x_mq <- read_dfdr_maxquant_msms(
  path = mq_path,
  pep_mode = "sanitize",          # or "drop" if PEP contains values >1
  exclude_contaminants = TRUE,
  add_decoy = 1L,
  unit = "psm",
  scope = "global",
  provenance = list(tool = "MaxQuant", file = basename(mq_path))
)
```

Next, all diagnostics measures can be obtained with:
```
diag <- dfdr_run_all(
  xs = list(MaxQuant_PSM = x_mq),
  alpha_main = 0.01,
  alphas = c(1e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1),
  eps = 0.2,
  win_rel = 0.2,
  truncation = "warn_drop",
  low_conf = c(0.2, 0.5),
  compute_pseudo_pvalues = TRUE  # <-- This adds p-value diagnostics
)
```

Finally, a single HTML report can be obtained using:
```
dfdr_render_report(diag, out_dir = "path/to/diagFDR_diann_out")
```

---

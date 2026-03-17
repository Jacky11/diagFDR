<h1 align="center">diagFDR: Verifiable FDR Diagnostics for Proteomics</h1>


<div align="center" style="margin-top: 10px;">
  <a href="https://www.gnu.org/licenses/gpl-3.0"><img src="https://img.shields.io/badge/License-GPLv3-blue.svg" alt="GPL v3 License"></a>
</div></a>
  
## 1. Description

False discovery rate (FDR) control is central to the credibility of peptide and protein identifications in mass spectrometry–based proteomics. **diagFDR** provides pipeline-agnostic diagnostics for target–decoy workflows that assess calibration and stability of reported confidence measures. It checks the coherence of scores, q-values, and posterior error probabilities (PEPs), quantifies decoy support near operating thresholds, measures cutoff sensitivity, and evaluates the equal-chance assumption using q-value band diagnostics.

## 2. Installation

It is recommended to install the latest version of R. 

The installation of the development version of the `diagFDR` package can be performed by:
```
if(!requireNamespace("devtools", quietly = TRUE)){
   install.packages("devtools") 
}
devtools::install_github("Jacky11/diagFDR")
```
You can now load it inside the R session:
```
library(diagFDR)
```

## 3. Verifiable FDR diagnostics HTML reports

A standardized HTML report can be obtained from any identification software using the dfdr_render_report() function to check the scope, the calibration, and the stability of a FDR-based identification search. Here we provide the way to get it from most common software, and from any software output.

### 3.1 DIA-NN

To enable all diagnostics at the precursor level, perform the search in DIA-NN using:

- decoys export: `--report-decoys`
- a permissive export ceiling: `--qvalue 0.5` (or higher if possible)

The q-value ceiling matters because some diagnostics operate in low-confidence regions
(e.g. equal-chance plausibility checks, or local-window support around cutoffs). 
The DIA-NN search has to provide you a `.parquet ` file containing the columns `Precursor.Id`, `Decoy`, `Q.Value` (optional: `Run`, `Global.Q.Value`, `PEP`).

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
  xs = list(global = x_global_gq, runx = x_runx),
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

Alternatively, the paper associated to the package proposes a common scope misuse (min-run aggregation). This can be performed by defining the universe with:
```
# Scope misuse comparator: min run-wise q over runs per precursor (anti-pattern)
x_minrun <- diann_global_minrunq(
  rep,
  q_col = "Q.Value",
  q_max_export = 0.5,
  unit = "precursor",
  scope = "aggregated",
  q_source = "min_run(Q.Value)"
)
```
and next by computing the diagnostic measures with `dfdr_run_all` as shown before.

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
dfdr_render_report(diag, out_dir = "path/to/diagFDR_maxquant_out")
```

### 3.3 mzIdentML files (Comet, X!Tandem, OMSSA, MS-GF+, Mascot, PEAKS, PeptideShaker)

Many search engines can export identifications in mzIdentML format. When explicit q-values
and/or PEPs are not present in the mzIdentML, `diagFDR` can:

1. extract a **competed PSM universe** (rank-1 by default),
2. infer target/decoy labels,
3. select a primary numeric score CV term (configurable),
4. reconstruct **TDC q-values** from scores via target--decoy counting (TDC),
5. compute scope/calibration/stability diagnostics.

mzIdentML is flexible and tool-dependent. If the score CV term or decoy labeling is not encoded consistently, you may need to adjust `score_accession_preference`, `score_direction`, and/or `decoy_regex`.

First, you have to upload the file in your R session:
```
mzid_path <- "path/to/search_results.mzid"

#Read the mzid file. Check the help of the function using `help(read_dfdr_mzid)` to adapt to your software outputs. 
x_mzid <- read_dfdr_mzid(
  mzid_path = mzid_path,
  rank = 1L,  # competed universe: take rank-1 SpectrumIdentificationItem

  # Choose a score CV term (priority list) and interpret its direction
  score_accession_preference = c(
    "MS:1002257", # example: MS-GF:RawScore (often higher is better)
    "MS:1001330", # Mascot:score (higher is better)
    "MS:1001328", # SEQUEST:xcorr (higher is better)
    "MS:1002052",
    "MS:1002049",
    "MS:1001331",
    "MS:1001171",
    "MS:1001950",
    "MS:1002466"
  ),
  score_direction = "auto",  # or "higher_better"/"lower_better" if auto fails

  # TDC correction: FDR_hat = (D + add_decoy)/T
  add_decoy = 1L,

  # Strict by default: require score for all PSMs (set <1 to allow missing)
  min_score_coverage = 1.0,

  # Fallback decoy inference if PeptideEvidence@isDecoy is not informative
  decoy_regex = "(^##|_REVERSED$|^REV_|^DECOY_)",

  unit = "psm",
  scope = "global",
  provenance = list(file = basename(mzid_path))
)
```

Next, all diagnostics measures can be obtained with:
```
diag <- dfdr_run_all(
  xs = list(mzid_PSM = x_mzid),
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
dfdr_render_report(diag, out_dir = "path/to/diagFDR_mzid_out")
```
### 3.4 Spectronaut

The code below shows how to run **diagFDR** directly from an export of Spectronaut at the "elution group" level. 
The export is assumed to be a "normal" report (not "pivot" report) that contains peptides characterized by "elution groups". 
It has to contain the columns "EG.Qvalue", "EG.Cscore" and "EG.IsDecoy".

First, you have to upload the file in your R session:
```
rep <- read_spectronaut_efficient("path/to/search_results.Report-Peptide normal.tsv", minimal = TRUE, dec = ",")
```
Next, the elution group universe can be extracted by:
```
univ_runwise <- spectronaut_runxprecursor(
  rep,
  q_col = "EG.Qvalue",
  score_col = "EG.Cscore"
)
```
Next, all diagnostics measures can be obtained with:
```
diag <- dfdr_run_all(
  xs = list(runwise = univ_runwise),
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
dfdr_render_report(diag, out_dir = "path/to/diagFDR_spectronaut_out")
```

### 3.5 Any software output

Any software producing outputs that can be mapped to columns `id`, `is_decoy`, `q`, `pep`, `run`, and `score` can similarly be used.
Here, we provide a simulated dataset:
```
n <- 8000
toy <- data.frame(
  id = as.character(seq_len(n)),
  is_decoy = sample(c(FALSE, TRUE), n, replace = TRUE, prob = c(0.98, 0.02)),
  run = sample(paste0("run", 1:3), n, replace = TRUE),
  score = c(rnorm(n * 0.98, mean = 7, sd = 1), rnorm(n * 0.02, mean = 5, sd = 1))[seq_len(n)],
  pep = NA_real_
)

# Reconstruct q-values by simple TDC from score (for toy data only):
# Here we mimic a typical "higher score is better" setting.
toy <- toy[order(toy$score, decreasing = TRUE), ]
toy$D_cum <- cumsum(toy$is_decoy)
toy$T_cum <- cumsum(!toy$is_decoy)
toy$FDR_hat <- (toy$D_cum + 1) / pmax(toy$T_cum, 1)
toy$q <- rev(cummin(rev(toy$FDR_hat)))
toy <- toy[, c("id","is_decoy","q","pep","run","score")]
```

Next, the dataset has to be converted into a `dfdr_tbl`:
```
x_toy <- as_dfdr_tbl(
  toy,
  unit = "psm",
  scope = "global",
  q_source = "toy TDC from score",
  q_max_export = 1,
  provenance = list(tool = "toy")
)
```
Next, all diagnostics measures can be obtained with:
```
diag <- dfdr_run_all(
  xs = list(univ = x_toy),
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
dfdr_render_report(diag, out_dir = "path/to/diagFDR_sim_out")
```

---

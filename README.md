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

## 3. 

---

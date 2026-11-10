# Example dataset

A simulated label-free experiment in the *Generic* wide format, read by the "Load example dataset"
button of the app and used in the vignettes: 1,500 proteins and 800 phosphosites, two conditions
(control, treated) with four samples each, the same samples in both tables.

* `protein.tsv`: `protein_id`, `gene`, `intensity_<sample>` (raw intensities, NA = missing)
* `ptm.tsv`: `feature_id`, `protein_ids`, `gene`, `stripped_sequence`, `modified_sequence`,
  `residue`, `position`, `localization_prob`, `ptm_type`, `intensity_<sample>`
* `annotation.tsv`: `sample`, `condition`

Generated with `simulate_ptm_experiment()` (seed 1); 10 % of sites change in usage, 10 % with
their protein only, 5 % are compensated (usage changes, intensity does not).

toy_diann <- function(n = 6000) {
  set.seed(2)
  dec <- rep(c(TRUE, FALSE), c(n / 3, 2 * n / 3))
  # decoy-rich low-confidence band to create an equal-chance imbalance
  q <- c(runif(n / 3, 0.2, 0.5), c(runif(n / 3, 0, 0.01), runif(n / 3, 0.2, 0.5)))
  data.frame(Precursor.Id = paste0("p", seq_len(n)), Run = "r1",
             Decoy = as.integer(dec), Q.Value = q, Global.Q.Value = q,
             PEP = pmin(1, q * 2), stringsAsFactors = FALSE)
}

test_that("a DIA-NN list without two_pass gives a caution, not an alert", {
  rep <- toy_diann()
  x <- suppressMessages(diann_global_precursor(rep, q_col = "Global.Q.Value", q_max_export = 0.5))
  h <- flag_headline(data.frame(list = "a", alpha = 0.01, effect_abs = 0.2,
                                eq_evaluable = TRUE, eq_reason = "pass_unknown"))
  expect_equal(h$flag_equalchance, "CAUTION_imbalance_pass_unknown")
  expect_equal(h$status, "CAUTION")
  h2 <- flag_headline(data.frame(list = "a", alpha = 0.01, effect_abs = 0.2,
                                 eq_evaluable = FALSE, eq_reason = "two_pass"))
  expect_equal(h2$flag_equalchance, "NE_two_pass")
  h3 <- flag_headline(data.frame(list = "a", alpha = 0.01, effect_abs = 0.2,
                                 eq_evaluable = TRUE, eq_reason = NA_character_))
  expect_equal(h3$flag_equalchance, "ALERT_imbalance")
})

test_that("dfdr_run_all marks the pass status of DIA-NN lists", {
  rep <- toy_diann()
  x_na <- suppressMessages(diann_global_precursor(rep, q_col = "Global.Q.Value", q_max_export = 0.5))
  x_tp <- diann_global_precursor(rep, q_col = "Global.Q.Value", q_max_export = 0.5, two_pass = TRUE)
  d <- suppressWarnings(dfdr_run_all(list(na = x_na, tp = x_tp), alpha_main = 0.01,
                                     compute_pseudo_pvalues = FALSE, min_N_equalchance = 100))
  h <- d$tables$headline
  expect_equal(h$eq_reason[h$list == "na"], "pass_unknown")
  expect_equal(h$eq_reason[h$list == "tp"], "two_pass")
})

test_that("the report renders with a relative out_dir and has table legends", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available())
  rep <- toy_diann()
  x <- diann_global_precursor(rep, q_col = "Global.Q.Value", q_max_export = 0.5, two_pass = TRUE)
  d <- suppressWarnings(dfdr_run_all(list(g = x), alpha_main = 0.01, compute_pseudo_pvalues = FALSE))
  tmp <- withr::local_tempdir()
  withr::local_dir(tmp)
  f <- dfdr_render_report(d, out_dir = "rep_out")
  expect_true(file.exists(f))
  html <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_match(html, "table-legend", fixed = TRUE)
  expect_false(grepl("$D_\\alpha", html, fixed = TRUE))
})

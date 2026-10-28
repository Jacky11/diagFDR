test_that("dfdr_ipe uses the D_b/T_b comparator", {
  set.seed(2)
  n <- 40000
  null <- runif(n) < 0.4
  dec <- null & runif(n) < 0.5
  pep <- ifelse(null, runif(n, 0, 0.5), runif(n, 0, 0.05))
  x <- as_dfdr_tbl(tibble::tibble(id = as.character(seq_len(n)), is_decoy = dec,
                                  q = 0, pep = pep, run = "r", score = NA_real_),
                   unit = "psm", scope = "global", q_source = "toy")
  r <- dfdr_ipe(x, n_min = 50)
  b <- r$bins[r$bins$used, ]
  w <- b$n / sum(b$n)
  expect_equal(r$IPE, sum(w * abs(b$pep_mean_targets - pmin(1, b$n_decoy / b$n_target))))
  expect_true(r$bins_used > 0)
})

test_that("flags are alerts and equal chance is n/e on two-pass exports", {
  h <- tibble::tibble(list = c("a", "b"), alpha = 0.01, D_alpha = c(5, 800),
                      support_index = 1 / sqrt(c(5, 800)), D_alpha_win = c(1, 100),
                      IPE = c(0.075, 0.01), FDR_hat = c(0.0223, 0.01),
                      effect_abs = c(0.2, 0.2), eq_evaluable = c(TRUE, FALSE),
                      eq_reason = c(NA, "two_pass"))
  f <- flag_headline(h)
  expect_equal(f$flag_IPE, c("ALERT_above_0.05", "OK"))
  expect_equal(f$flag_equalchance, c("ALERT_imbalance", "NE_two_pass"))
  expect_equal(f$status, c("ALERT", "NO_ALERT"))
  expect_equal(f$flag_FDR[1], "ALERT_ratio_above_alpha")
})

test_that("headline reports support index and FDP-CV benchmark", {
  x <- as_dfdr_tbl(tibble::tibble(id = as.character(1:200),
                                  is_decoy = rep(c(FALSE, TRUE), c(150, 50)),
                                  q = 0.001, pep = NA_real_, run = "r", score = NA_real_),
                   unit = "psm", scope = "global", q_source = "toy")
  h <- dfdr_headline(x, alpha = 0.01)
  expect_equal(h$support_index, 1 / sqrt(50))
  expect_equal(h$CV_benchmark, sqrt(2 / 50))
  
})

test_that("mokapot recount gives (D+1)/T q-values on winners", {
  raw <- tibble::tibble(run = "r", spectrum_id = as.character(1:6),
                        `mokapot score` = c(6, 5, 4, 3, 2, 1),
                        `mokapot q-value` = 0.5, `mokapot PEP` = 0.1,
                        is_target = c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE))
  x <- mokapot_competed_universe(raw, q_mode = "recount")
  expect_equal(attr(x, "meta")$provenance$q_mode, "recount")
  expect_true(all(x$q <= 1))
  expect_equal(x$q[x$id == "r||1"], 0.5)
})

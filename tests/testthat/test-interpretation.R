toy_comp <- function(n = 20000, seed = 3) {
  set.seed(seed)
  null <- runif(n) < 0.6
  xt <- rnorm(n) + ifelse(null, 0, 3); yd <- rnorm(n)
  win <- pmax(xt, yd); dec <- yd > xt
  q <- diagFDR:::tdc_qvalues(win, dec, add_decoy = 1L)
  f0 <- dnorm(win); f1 <- dnorm(win - 3)
  pep <- pmin(1, 0.6 * f0 / (0.6 * f0 + 0.4 * f1))
  data.frame(id = paste0("u", seq_len(n)), is_decoy = dec, q = q, pep = pep,
             run = "r", score = win)
}

test_that("step-up threshold reproduces Benjamini-Hochberg", {
  set.seed(1)
  p <- c(runif(900), rbeta(100, 0.1, 5))
  t <- diagFDR:::.stepup_threshold(p, 0.05, 1)
  expect_equal(sum(p <= t), sum(stats::p.adjust(p, "BH") <= 0.05))
  expect_true(is.na(diagFDR:::.stepup_threshold(runif(100, 0.5, 1), 0.01, 1)))
})

test_that("dfdr_p_summary gives pi0 by subset and thresholds", {
  df <- toy_comp()
  df$p <- score_to_pvalue(df$score, method = "decoy_ecdf", is_decoy = df$is_decoy)
  x <- as_dfdr_tbl(df, unit = "psm", scope = "global", q_source = "toy")
  attr(x, "meta")$p_source <- "score_to_pvalue(score, method='decoy_ecdf')"
  s <- dfdr_p_summary(x, alpha = 0.01)
  expect_setequal(s$summary$subset, c("all", "target", "decoy"))
  th <- s$thresholds
  est <- th[th$subset == "target" & th$pi0_choice == "estimate", ]
  bh <- th[th$subset == "target" & th$pi0_choice == "1 (BH)", ]
  expect_gte(est$p_threshold, bh$p_threshold)
  expect_equal(est$n_accepted_decoy, 0L)
  expect_true(is.finite(est$score_threshold))
  txt <- dfdr_p_interpretation(s, T_alpha = sum(!df$is_decoy & df$q <= 0.01))
  expect_type(txt, "character")
  expect_true(any(grepl("by construction", txt)))
  expect_true(any(grepl("score", txt)))
})

test_that("too few p-values are reported as such", {
  df <- toy_comp(n = 500)
  df$p <- runif(500)
  x <- as_dfdr_tbl(df, unit = "psm", scope = "global", q_source = "toy")
  s <- dfdr_p_summary(x, alpha = 0.01, min_n = 2000)
  expect_match(dfdr_p_interpretation(s)[1], "Too few p-values")
})

test_that("PEP summary detects under- and overestimation", {
  df <- toy_comp()
  x_ok <- as_dfdr_tbl(df, unit = "psm", scope = "global", q_source = "toy")
  df_u <- df; df_u$pep <- df_u$pep / 4
  df_o <- df; df_o$pep <- pmin(1, df_o$pep * 4)
  su <- dfdr_pep_summary(as_dfdr_tbl(df_u, unit = "psm", scope = "global", q_source = "toy"))
  so <- dfdr_pep_summary(as_dfdr_tbl(df_o, unit = "psm", scope = "global", q_source = "toy"))
  expect_equal(su$verdict, "under")
  expect_equal(so$verdict, "over")
  expect_match(dfdr_pep_interpretation(su), "underestimate")
  expect_match(dfdr_pep_interpretation(so), "overestimate")
  expect_true(dfdr_pep_summary(x_ok)$verdict %in% c("agree", "under", "over", "mixed"))
})

test_that("pi0 consensus groups converging estimators and lists the others", {
  m <- c("st.spline", "st.boot", "jiang", "histo", "langaas", "pounds", "abh", "slim")
  v <- c(0.30, 0.31, 0.29, 0.31, 0.283, 0.34, 0.877, 0.62)
  cs <- diagFDR:::.pi0_consensus(v, m)
  expect_true(cs$converged)
  expect_equal(cs$n_in, 6L)
  expect_setequal(cs$others, c("abh", "slim"))
  expect_equal(cs$value, stats::median(v[1:6]))
  cs2 <- diagFDR:::.pi0_consensus(c(0.1, 0.4, 0.7, 0.9), m[1:4])
  expect_false(cs2$converged)
})

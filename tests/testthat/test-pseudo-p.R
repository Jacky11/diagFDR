mk <- function(with_score = TRUE, n = 400) {
  set.seed(1)
  dec <- rep(c(TRUE, FALSE), each = n / 2)
  sc <- c(rnorm(n / 2), rnorm(n / 2, 2))
  df <- data.frame(id = paste0("u", seq_len(n)), is_decoy = dec,
                   q = runif(n), pep = plogis(-sc), run = "r1",
                   score = if (with_score) sc else NA_real_)
  as_dfdr_tbl(df, unit = "precursor", scope = "global", q_source = "toy")
}

test_that("pseudo-p-values are computed from the score", {
  out <- diagFDR:::add_pseudo_pvalues_to_list(list(a = mk(TRUE)))
  expect_true(any(is.finite(out$a$p)))
  expect_match(attr(out$a, "meta")$p_source, "score_to_pvalue(score", fixed = TRUE)
})

test_that("a list without score warns and records why", {
  expect_warning(out <- diagFDR:::add_pseudo_pvalues_to_list(list(a = mk(FALSE))),
                 "no score column")
  expect_false("p" %in% names(out$a) && any(is.finite(out$a$p)))
  expect_match(attr(out$a, "meta")$p_source, "not computed")
})

test_that("the PEP fallback gives pseudo-p-values", {
  out <- diagFDR:::add_pseudo_pvalues_to_list(list(a = mk(FALSE)), score_source = "pep")
  expect_true(any(is.finite(out$a$p)))
  expect_match(attr(out$a, "meta")$p_source, "-PEP")
})

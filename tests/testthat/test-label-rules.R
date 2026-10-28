test_that("the four label rules reproduce the worked example of Note S9", {
  rep <- data.frame(Protein.Group = rep("G1", 4), Decoy = c(1, 0, 0, 0),
                    Global.PG.Q.Value = c(0.004, 0.02, 0.30, 0.60))
  u <- function(r) collapse_universe(rep, "Protein.Group", "Global.PG.Q.Value", rule = r)
  a <- u("any_decoy");  expect_equal(nrow(a), 1); expect_true(a$is_decoy);  expect_equal(a$q, 0.004)
  w <- u("winning_row"); expect_equal(nrow(w), 1); expect_true(w$is_decoy)
  m <- u("majority");   expect_equal(nrow(m), 1); expect_false(m$is_decoy); expect_equal(m$q, 0.004)
  s <- u("split");      expect_equal(nrow(s), 2)
  expect_equal(sort(s$q), c(0.004, 0.02))
  expect_equal(s$is_decoy[s$q == 0.004], TRUE)
  expect_equal(attr(s, "meta")$provenance$n_mixed, 1)
})

test_that("rules coincide when the label is intrinsic", {
  rep <- data.frame(Precursor.Id = c("A", "A", "B", "C"), Decoy = c(0, 0, 1, 0),
                    Q.Value = c(0.01, 0.002, 0.3, 0.05))
  expect_true(label_is_intrinsic(rep, "Precursor.Id")$intrinsic)
  out <- dfdr_label_rule_audit(rep, "Precursor.Id", "Q.Value", alpha = 0.05)
  expect_equal(length(unique(out$D_alpha)), 1)
  expect_equal(length(unique(out$T_alpha)), 1)
})

test_that("factor Decoy columns are not misread as level codes", {
  rep <- data.frame(Protein.Group = c("G1", "G2"), Decoy = factor(c("1", "0")),
                    Global.PG.Q.Value = c(0.1, 0.2))
  s <- collapse_universe(rep, "Protein.Group", "Global.PG.Q.Value")
  expect_equal(s$is_decoy[s$id == "G1"], TRUE)
  expect_equal(s$is_decoy[s$id == "G2"], FALSE)
})

test_that("diann_global_proteingroup wraps collapse_universe", {
  rep <- data.frame(Protein.Group = c("P1", "P1", "P2"), Decoy = c(1, 0, 0),
                    Global.PG.Q.Value = c(0.2, 0.001, 0.003))
  x <- diann_global_proteingroup(rep, two_pass = FALSE)
  expect_s3_class(x, "dfdr_tbl")
  expect_equal(nrow(x), 3)
  expect_equal(attr(x, "meta")$unit, "protein_group")
})

test_that("raw and corrected entrapment estimates reproduce the paper's arithmetic", {
  ## Section 4.2: T = 81,643, E = 177, r = 0.432 -> raw 0.217%, corrected 0.719%
  x <- as_dfdr_tbl(data.frame(id = as.character(1:81643), is_decoy = FALSE,
                              q = 0.001, pep = NA_real_, run = NA_character_, score = NA_real_))
  ent <- as.character(1:177)
  e <- dfdr_entrapment(x, ent, alpha = 0.01, r = 0.432)
  expect_equal(e$E_alpha, 177); expect_equal(e$T_alpha, 81643)
  expect_equal(round(100 * e$fdp_raw, 3), 0.217)
  expect_equal(round(100 * e$fdp_corrected, 3), 0.719)
  expect_equal(round(100 * e$fdp_corrected_target_only, 3), 0.503)
  expect_equal(round(100 * e$lo, 3), 0.187); expect_equal(round(100 * e$hi, 3), 0.251)
})

test_that("discordant-set rate separates core and added sets", {
  a1 <- paste0("p", 1:100); a2 <- paste0("p", 1:120)
  d <- dfdr_entrapment_discordant(a1, a2, entrapment_ids = c("p3", "p101", "p115"))
  expect_equal(d$E[d$set == "core"], 1)
  expect_equal(d$n[d$set == "added"], 20); expect_equal(d$E[d$set == "added"], 2)
  expect_equal(d$n[d$set == "removed"], 0)
  expect_true(attr(d, "nested"))
  expect_equal(attr(d, "ratio"), (2 / 20) / (1 / 100))
})

test_that("entrapment attribution uses accession sets, not prefixes", {
  expect_equal(dfdr_is_entrapment_group(c("A1;A2", "A1;P1", "P1", NA), c("A1", "A2")),
               c(TRUE, FALSE, FALSE, FALSE))
  expect_equal(dfdr_is_entrapment_group(c("A1;A2", "A1;P1"), c("A1", "A2"), rule = "any"),
               c(TRUE, TRUE))
  rep <- data.frame(Precursor.Id = c("x", "y"), Protein.Group = c("P1", "P2"))
  expect_error(diann_entrapment_ids(rep, accessions = "ARATH"), "No identifier")
})

test_that("FASTA parsing and digestion give the expected size ratios", {
  t <- tempfile(); e <- tempfile()
  writeLines(c(">sp|P1|A", "MKAAAAAAAKRPEPTIDEKLLLLLLLR", ">sp|P2|B", "GGGGGGGGKRRRR"), t)
  writeLines(c(">sp|A1|C", "MKAAAAAAAKRQQQQQQQQK"), e)
  expect_equal(read_fasta_accessions(t), c("P1", "P2"))
  r <- dfdr_entrapment_size_ratio(t, e, missed = 0)
  expect_equal(r$value[r$quantity == "r_protein"], 0.5)
  expect_equal(r$value[r$quantity == "shared_peptides"], 1)   # AAAAAAAK
})

test_that("picked-protein competition keeps the better of each pair and recounts", {
  df <- data.frame(id = c("P1", "P1", "P2", "P3", "P3"),
                   is_decoy = c(FALSE, TRUE, FALSE, TRUE, FALSE),
                   q = c(0.001, 0.2, 0.002, 0.004, 0.05))
  p <- dfdr_picked_protein(df, alpha = 0.5)
  expect_equal(p$n_pairs, 2)
  expect_equal(nrow(p$tbl), 3)
  expect_equal(p$tbl$is_decoy[p$tbl$id == "P3"], TRUE)   # decoy 0.004 beats target 0.05
  expect_equal(p$n_pick_T, 2); expect_equal(p$n_pick_D, 1)
  expect_true(all(diff(p$tbl$q[order(p$tbl$q_in)]) >= 0))   # monotone q
})

test_that("picked competition agrees with split on a list without mixed groups", {
  sim <- dfdr_simulate_tdc(n = 5000, seed = 11)
  p <- dfdr_picked_protein(sim, alpha = 0.01)
  h <- dfdr_headline(sim, alpha = 0.01)
  expect_equal(p$T_alpha, h$T_alpha)
  expect_equal(p$D_alpha, h$D_alpha)
})

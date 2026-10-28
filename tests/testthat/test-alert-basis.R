test_that("dfdr_alert_basis lists every flag with a basis", {
  b <- dfdr_alert_basis()
  expect_true(all(c("flag_FDR", "flag_Dalpha", "flag_support", "flag_Dwin",
                    "flag_IPE", "flag_equalchance") %in% b$flag))
  expect_true(all(b$basis %in% c("derived check", "simulation-informed alert",
                                 "heuristic alert", "descriptive")))
  expect_equal(b$basis[b$flag == "flag_IPE"], "heuristic alert")
})

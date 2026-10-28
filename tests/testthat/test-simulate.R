test_that("the simulator holds equal chance and the decoy ratio at alpha", {
  sim <- dfdr_simulate_tdc(n = 50000, seed = 5)
  f <- dfdr_realised_fdp(sim, alpha = 0.01)
  expect_true(abs(f$decoy_ratio - 0.01) < 0.002)
  expect_true(f$FDP > 0 && f$FDP < 0.03)
  low <- sim$q >= 0.2 & sim$q <= 0.5
  expect_true(abs(mean(sim$is_decoy[low]) - 0.5) < 0.05)
})

test_that("decoy misspecification inflates the realised FDP but not the decoy ratio", {
  s1 <- dfdr_simulate_tdc(n = 50000, decoy_scale = 0.85, seed = 6)
  f <- dfdr_realised_fdp(s1, alpha = 0.01)
  expect_true(abs(f$decoy_ratio - 0.01) < 0.002)
  expect_true(f$FDP > 0.015)
})

test_that("kappa is near sqrt(2) under independent Gaussian nulls", {
  skip_on_cran()
  k <- dfdr_check_cv_benchmark(B = 40, n = 20000, alpha = 0.05)
  expect_true(k$kappa > 1.1 && k$kappa < 1.8)
})

test_that("the bundled example runs through the headline", {
  h <- dfdr_headline(sim_tdc_example, alpha = 0.01)
  expect_true(h$D_alpha > 0)
})

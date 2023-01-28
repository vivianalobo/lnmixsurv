test_that("print method works", {
  mod <- new_survival_ln_fit_with_fixed_rng()
  expect_snapshot(mod)
})

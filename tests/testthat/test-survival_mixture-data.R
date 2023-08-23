mod_spec <- parsnip::survival_reg(mode = "censored regression", engine = "survival_ln_mixture")
f_fit <- withr::with_rng_version(
  "4.2.2",
  withr::with_seed(
    1,
    parsnip::fit(mod_spec, survival::Surv(y, delta) ~ x, data = sim_data$data) # nolint: object_usage_linter.
  )
)

mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))

test_that("parsnip specification works", {
  expect_equal(f_fit$fit, mod, tolerance = 10^-1)
})

test_that("parsnip survival prediction works", {
  new_data <- data.frame(x = c("0", "1"))
  pred <- predict(mod, new_data = new_data, type = "survival", eval_time = c(20, 100), interval = "credible", level = 0.8)
  expected <- predict(f_fit, new_data = new_data, type = "survival", eval_time = c(20, 100), interval = "credible", level = 0.8)

  expect_equal(pred, expected, tolerance = 10^-1)
})

test_that("parsnip hazard prediction works", {
  new_data <- data.frame(x = c("0", "1"))
  pred <- predict(mod, new_data = new_data, type = "hazard", eval_time = c(20, 100))
  expected <- predict(f_fit, new_data = new_data, type = "hazard", eval_time = c(20, 100))

  expect_equal(pred, expected, tolerance = 10^-1)
})

test_that("parsnip wont allow hazard predictions to have a interval", {
  new_data <- data.frame(x = c("0", "1"))
  pred <- predict(mod, new_data = new_data, type = "hazard", eval_time = c(20, 100), interval = "none")
  expected <- predict(f_fit, new_data = new_data, type = "hazard", eval_time = c(20, 100), interval = "credible", level = 0.8)
  
  expect_equal(pred, expected, tolerance = 10^-1)
})

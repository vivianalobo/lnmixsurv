test_that("outcome must be Surv object", {
  expect_error(
    survival_ln_mixture(y ~ x, sim_data$data)
  )
})

test_that("only right censored data supported", {
  expect_error(
    survival_ln_mixture(survival::Surv(y, delta, type = "left") ~ x, sim_data$data)
  )
})

test_that("events at time zero are ano supported", {
  data_new <- sim_data$data
  data_new$y[1] <- 0
  expect_error(
    survival_ln_mixture(survival::Surv(y, delta, type = "left") ~ x, sim_data$data))
})

test_that("when using ~NULL, intercept must be TRUE", {
  expect_error(
    survival_ln_mixture(survival::Surv(y, delta) ~ NULL, sim_data$data, intercept = FALSE)
  )
})

test_that("survival_ln_mixture doesnt work with xy specification", {
  expect_error(
    survival_ln_mixture(sim_data$data$y, sim_data$data$x)
  )
})

test_that("survival_ln_mixture works with intercept only fit", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_intercept_only.rds"))
  expect_equal(tidy(mod)$estimate, c(3.882219, 4.875943), 
               tolerance = 10^-1)
})

test_that("fit works as expected with simulated data", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  post_summary <- posterior::summarise_draws(mod$posterior, estimate = stats::median, std.error = stats::mad)
  colnames(post_summary)[1] <- "term"
  post_tidy <- tidy(mod, effects = c("fixed", "auxiliary"))
  expected_result <- structure(
    list(term = c(
      "(Intercept)_1", "x1_1", "(Intercept)_2",
      "x1_2", "phi_1", "phi_2", "eta_1"
    ), estimate = c(
      4.0449468,
      0.8091331, 3.4289818, 0.4880462, 26.8527903,
      3.1718785, 0.5048220
    ), std.error = c(
      0.006640402,
      0.009723673, 0.017118112, 0.021268204, 1.250917196,
      0.090543019, 0.009105166
    )),
    row.names = c(NA, -7L), class = c("draws_summary", "tbl_df", "tbl", "data.frame"), num_args = list()
  )

  expect_equal(mod$nobs, 10000)
  expect_equal(post_summary, post_tidy, tolerance = 10^-1)
  expect_equal(post_summary, expected_result, tolerance = 10^-1)
})

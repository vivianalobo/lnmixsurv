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
  expect_equal(tidy(mod)$estimate, c(4.345570, 3.434564),
               tolerance = 10^-2)
})

test_that("fit works as expected with simulated data", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  post_summary <- posterior::summarise_draws(mod$posterior, estimate = stats::median, std.error = stats::mad)
  colnames(post_summary)[1] <- "term"
  post_tidy <- tidy(mod, effects = c("fixed", "auxiliary"))
  expected_result <- structure(
    list(term = c(
      "(Intercept)_a", "x1_a", "(Intercept)_b",
      "x1_b", "phi_a", "phi_b", "theta_a"
    ), estimate = c(
      4.0445757,
      0.8090734, 3.4245488, 0.4915648, 26.3862345,
      3.1816509, 0.5088579
    ), std.error = c(
      0.006236653,
      0.009518335, 0.019522311, 0.021621560, 1.34253028,
      0.10961297, 0.01305763
    )),
    row.names = c(NA, -7L), class = c("draws_summary", "tbl_df", "tbl", "data.frame"), num_args = list()
  )

  expect_equal(mod$nobs, 10000, tolerance = 10^-2)
  expect_equal(post_summary, post_tidy, tolerance = 10^-2)
  expect_equal(post_summary, expected_result, tolerance = 10^-2)
})

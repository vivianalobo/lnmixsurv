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
  expect_equal(tidy(mod)$estimate, c(3.667312, 4.521098),
               tolerance = 10^-2)
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
      4.0446964,
      0.8095789, 3.4277493, 0.4868316, 26.5640273,
      3.1826322, 0.5051077
    ), std.error = c(
      0.006796203,
      0.009936305, 0.018647274, 0.021226930, 1.210941881,
      0.114646819, 0.011923455
    )),
    row.names = c(NA, -7L), class = c("draws_summary", "tbl_df", "tbl", "data.frame"), num_args = list()
  )

  expect_equal(mod$nobs, 10000)
  expect_equal(post_summary, post_tidy, tolerance = 10^-2)
  expect_equal(post_summary, expected_result, tolerance = 10^-2)
})

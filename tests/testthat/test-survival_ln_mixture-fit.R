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
  expect_equal(tidy(mod)$estimate, c(1.457353, 1.233633))
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
      1.3956390,
      0.1806862, 1.2153782, 0.1479326, 487.8716967,
      35.5621535, 0.5103092
    ), std.error = c(
      0.001518476,
      0.002307636, 0.004641445, 0.006147595, 21.18129631,
      1.04788310, 0.01014782
    )),
    row.names = c(NA, -7L), class = c("draws_summary", "tbl_df", "tbl", "data.frame"), num_args = list()
  )

  expect_equal(mod$nobs, 10000)
  expect_equal(post_summary, post_tidy)
  expect_equal(post_summary, expected_result)
})

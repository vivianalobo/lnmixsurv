test_that("tiders work as expected", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  obtained <- tidy(mod, conf.int = TRUE)
  expected <- structure(
    list(
      term = c("(Intercept)_a", "x1_a", "(Intercept)_b", "x1_b"),
      estimate = c(4.0445757, 0.8090734, 3.4245488, 0.4915648),
      std.error = c(0.006236653, 0.009518335, 0.019522311, 0.021621560),
      conf.low = c(4.0368830, 0.7958502, 3.4006760, 0.4627559),
      conf.high = c(4.0526401, 0.8199592, 3.4501877, 0.5175077)
    ), 
    row.names = c(NA, -4L), class = c("draws_summary", "tbl_df", "tbl", "data.frame"),
    num_args = list()
  )
  expect_equal(obtained, expected, tolerance = 10^-3)
})

test_that("tiders work as expected", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  obtained <- tidy(mod, conf.int = TRUE)
  expected <- structure(
    list(
      term = c("(Intercept)_1", "x1_1", "(Intercept)_2", "x1_2"),
      estimate = c(4.0455659, 0.8072690, 3.4262658, 0.4874025),
      std.error = c(0.006567584, 0.009969275, 0.020706785, 0.023423408),
      conf.low = c(4.0368490, 0.7946591, 3.3997479, 0.4582146),
      conf.high = c(4.0547803, 0.8201348, 3.4528735, 0.5154283)
    ),
    row.names = c(NA, -4L), class = c("draws_summary", "tbl_df", "tbl", "data.frame"),
    num_args = list()
  )
  expect_equal(obtained, expected, tolerance = 10^-1)
})

test_that("tiders work as expected", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  obtained <- tidy(mod, conf.int = TRUE)
  expected <- structure(
    list(
      term = c("(Intercept)_1", "x1_1", "(Intercept)_2", "x1_2"),
      estimate = c(4.0446964, 0.8095789, 3.4277493, 0.4868316),
      std.error = c(0.006796203, 0.009936305, 0.018647274, 0.021226930),
      conf.low = c(4.0359851, 0.7969076, 3.4033308, 0.4593627),
      conf.high = c(4.0533994, 0.8222505, 3.4498217, 0.5118913)
    ), 
    row.names = c(NA, -4L), class = c("draws_summary", "tbl_df", "tbl", "data.frame"),
    num_args = list()
  )
  expect_equal(obtained, expected, tolerance = 10^-2)
})

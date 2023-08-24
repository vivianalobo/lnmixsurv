test_that("tiders work as expected", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  obtained <- tidy(mod, conf.int = TRUE)
  expected <- structure(
    list(
      term = c("(Intercept)_1", "x1_1", "(Intercept)_2", "x1_2"),
      estimate = c(4.0456457, 0.8075922, 3.4202430, 0.4879574),
      std.error = c(0.006370159, 0.009649999, 0.020753769, 0.022196693),
      conf.low = c(4.0369899, 0.7948419, 3.3928802, 0.4613753),
      conf.high = c(4.0529927, 0.8188464, 3.4463243, 0.5157443)
    ), 
    row.names = c(NA, -4L), class = c("draws_summary", "tbl_df", "tbl", "data.frame"),
    num_args = list()
  )
  expect_equal(obtained, expected, tolerance = 10^-4)
})

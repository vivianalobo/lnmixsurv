test_that("tiders work as expected", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  obtained <- tidy(mod, conf.int = TRUE)
  expected <- structure(
    list(
      term = c("(Intercept)_a", "x1_a", "(Intercept)_b", "x1_b"),
      estimate = c(4.04383317528908, 0.808460445960304, 3.42274214257205, 0.490255319956356),
      std.error = c(0.00692710161402975, 0.00981968616603306, 0.021024321470229, 0.0193884470889992),
      conf.low = c(4.03558444537862, 0.796353441934959, 3.39210340879772, 0.461415419961591),
      conf.high = c(4.05375431039575, 0.821560789495629, 3.44910689813837, 0.513038499614571)
    ), 
    row.names = c(NA, -4L), class = c("draws_summary", "tbl_df", "tbl", "data.frame"),
    num_args = list()
  )
  expect_equal(obtained, expected)
})

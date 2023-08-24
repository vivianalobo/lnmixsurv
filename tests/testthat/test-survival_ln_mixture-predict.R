test_that("survival prediction works", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  new_data <- data.frame(x = c("0", "1"))
  pred <- predict(mod, new_data, type = "survival", eval_time = c(20, 100), interval = "credible")

  expected <- structure(
    list(.pred = list(structure(list(.eval_time = c(
      20,
      100
    ), .pred_survival = c(0.89078442, 0.00974128), 
    .pred_lower = c(0.885100441, 0.007757679), 
    .pred_upper = c(
      0.89646205,
      0.01220202
    )), row.names = c(NA, -2L), class = c(
      "tbl_df",
      "tbl", "data.frame"
    )), structure(list(
      .eval_time = c(20, 100),
      .pred_survival = c(0.9752529, 0.5107449),
      .pred_lower = c(0.9723125, 0.4994297), .pred_upper = c(
        0.9777107,
        0.5218468
      )
    ), row.names = c(NA, -2L), class = c(
      "tbl_df",
      "tbl", "data.frame"
    )))),
    class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -2L)
  )

  expect_equal(pred, expected, tolerance = 10^-1)
})

test_that("hazard prediction works", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  new_data <- data.frame(x = c("0", "1"))
  pred <- predict(mod, new_data = new_data, type = "hazard", eval_time = c(20, 100), interval = "credible")

  expected <- structure(list(.pred = list(structure(list(.eval_time = c(
    20,
    100
  ), .pred_hazard = c(0.01464948, 0.05653027), 
  .pred_lower = c(0.01405532, 0.04934615), 
  .pred_upper = c(
    0.01523202,
    0.06480105
  )), row.names = c(NA, -2L), class = c(
    "tbl_df",
    "tbl", "data.frame"
  )), structure(list(
    .eval_time = c(20, 100),
    .pred_hazard = c(0.004655651, 0.012265091),
    .pred_lower = c(0.004307983, 0.011422636),
    .pred_upper = c(0.005084482, 0.013103554)
  ), row.names = c(
    NA,
    -2L
  ), class = c("tbl_df", "tbl", "data.frame")))), class = c(
    "tbl_df",
    "tbl", "data.frame"
  ), row.names = c(NA, -2L))

  expect_equal(pred, expected, tolerance = 10^-1)
})

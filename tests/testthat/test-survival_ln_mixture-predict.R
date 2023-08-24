test_that("survival prediction works", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  new_data <- data.frame(x = c("0", "1"))
  pred <- predict(mod, new_data, type = "survival", eval_time = c(20, 100), interval = "credible")

  expected <- structure(
    list(.pred = list(structure(list(.eval_time = c(
      20,
      100
    ), .pred_survival = c(0.890715018, 0.009245273), 
    .pred_lower = c(0.88490914, 0.00722155), 
    .pred_upper = c(
      0.89675283,
      0.01174575
    )), row.names = c(NA, -2L), class = c(
      "tbl_df",
      "tbl", "data.frame"
    )), structure(list(
      .eval_time = c(20, 100),
      .pred_survival = c(0.9751677, 0.5110331),
      .pred_lower = c(0.9721102, 0.4984439), .pred_upper = c(
        0.9777439,
        0.5229762
      )
    ), row.names = c(NA, -2L), class = c(
      "tbl_df",
      "tbl", "data.frame"
    )))),
    class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -2L)
  )

  expect_equal(pred, expected, tolerance = 10^-4)
})

test_that("hazard prediction works", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  new_data <- data.frame(x = c("0", "1"))
  pred <- predict(mod, new_data = new_data, type = "hazard", eval_time = c(20, 100), interval = "credible")

  expected <- structure(list(.pred = list(structure(list(.eval_time = c(
    20,
    100
  ), .pred_hazard = c(0.01470163, 0.05774842), 
  .pred_lower = c(0.01406162, 0.05050136), 
  .pred_upper = c(
    0.01529979,
    0.06885545
  )), row.names = c(NA, -2L), class = c(
    "tbl_df",
    "tbl", "data.frame"
  )), structure(list(
    .eval_time = c(20, 100),
    .pred_hazard = c(0.004701174, 0.012289649),
    .pred_lower = c(0.00430829, 0.01139716),
    .pred_upper = c(0.005121054, 0.013121967)
  ), row.names = c(
    NA,
    -2L
  ), class = c("tbl_df", "tbl", "data.frame")))), class = c(
    "tbl_df",
    "tbl", "data.frame"
  ), row.names = c(NA, -2L))

  expect_equal(pred, expected, tolerance = 10^-4)
})

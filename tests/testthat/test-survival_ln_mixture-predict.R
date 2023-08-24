test_that("survival prediction works", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  new_data <- data.frame(x = c("0", "1"))
  pred <- predict(mod, new_data, type = "survival", eval_time = c(20, 100), interval = "credible")

  expected <- structure(
    list(.pred = list(structure(list(.eval_time = c(
      20,
      100
    ), .pred_survival = c(0.890899949, 0.009717829), 
    .pred_lower = c(0.885061278, 0.007802484), 
    .pred_upper = c(
      0.89698707,
      0.01228781
    )), row.names = c(NA, -2L), class = c(
      "tbl_df",
      "tbl", "data.frame"
    )), structure(list(
      .eval_time = c(20, 100),
      .pred_survival = c(0.9749848, 0.5093512),
      .pred_lower = c(0.9721684, 0.4983275), .pred_upper = c(
        0.9777895,
        0.5200507
      )
    ), row.names = c(NA, -2L), class = c(
      "tbl_df",
      "tbl", "data.frame"
    )))),
    class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -2L)
  )

  expect_equal(pred, expected)
})

test_that("hazard prediction works", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  new_data <- data.frame(x = c("0", "1"))
  pred <- predict(mod, new_data = new_data, type = "hazard", eval_time = c(20, 100), interval = "credible")

  expected <- structure(list(.pred = list(structure(list(.eval_time = c(
    20,
    100
  ), .pred_hazard = c(0.01465875, 0.05620879), 
  .pred_lower = c(0.01406332, 0.04897837), 
  .pred_upper = c(
    0.01526089,
    0.06593787
  )), row.names = c(NA, -2L), class = c(
    "tbl_df",
    "tbl", "data.frame"
  )), structure(list(
    .eval_time = c(20, 100),
    .pred_hazard = c(0.004698114, 0.012236965),
    .pred_lower = c(0.004305749, 0.011391576),
    .pred_upper = c(0.005108614, 0.013169467)
  ), row.names = c(
    NA,
    -2L
  ), class = c("tbl_df", "tbl", "data.frame")))), class = c(
    "tbl_df",
    "tbl", "data.frame"
  ), row.names = c(NA, -2L))

  expect_equal(pred, expected)
})

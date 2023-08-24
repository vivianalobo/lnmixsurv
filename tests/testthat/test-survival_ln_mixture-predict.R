test_that("survival prediction works", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  new_data <- data.frame(x = c("0", "1"))
  pred <- predict(mod, new_data, type = "survival", eval_time = c(20, 100), interval = "credible")

  expected <- structure(
    list(.pred = list(structure(list(.eval_time = c(
      20,
      100
    ), .pred_survival = c(0.891042225, 0.009844677), 
    .pred_lower = c(0.885530907, 0.007567918), 
    .pred_upper = c(
      0.89670611,
      0.01196521
    )), row.names = c(NA, -2L), class = c(
      "tbl_df",
      "tbl", "data.frame"
    )), structure(list(
      .eval_time = c(20, 100),
      .pred_survival = c(0.9750390, 0.5093482),
      .pred_lower = c(0.9719968, 0.4982551), .pred_upper = c(
        0.9774285,
        0.5204384
      )
    ), row.names = c(NA, -2L), class = c(
      "tbl_df",
      "tbl", "data.frame"
    )))),
    class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -2L)
  )

  expect_equal(pred, expected, tolerance = 10^-2)
})

test_that("hazard prediction works", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  new_data <- data.frame(x = c("0", "1"))
  pred <- predict(mod, new_data = new_data, type = "hazard", eval_time = c(20, 100), interval = "credible")

  expected <- structure(list(.pred = list(structure(list(.eval_time = c(
    20,
    100
  ), .pred_hazard = c(0.01469014, 0.05583335), 
  .pred_lower = c(0.01407298, 0.04952352), 
  .pred_upper = c(
    0.01528858,
    0.06514244
  )), row.names = c(NA, -2L), class = c(
    "tbl_df",
    "tbl", "data.frame"
  )), structure(list(
    .eval_time = c(20, 100),
    .pred_hazard = c(0.00470530, 0.01216557),
    .pred_lower = c(0.004359482, 0.011389595),
    .pred_upper = c(0.00512277, 0.01312007)
  ), row.names = c(
    NA,
    -2L
  ), class = c("tbl_df", "tbl", "data.frame")))), class = c(
    "tbl_df",
    "tbl", "data.frame"
  ), row.names = c(NA, -2L))

  expect_equal(pred, expected, tolerance = 10^-2)
})

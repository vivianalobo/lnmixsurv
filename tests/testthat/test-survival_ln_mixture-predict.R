test_that("survival prediction works", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  new_data <- data.frame(x = c("0", "1"))
  pred <- predict(mod, new_data, type = "survival", eval_time = c(20, 100), interval = "credible")

  expected <- structure(list(.pred = list(structure(list(
    .time = c(20, 100),
    .pred_survival = c(0.890693737049813, 0.00982157645355162), .pred_lower = c(0.885272530566163, 0.00772484808096847), .pred_upper = c(0.89664456777491, 0.0127486855311491)
  ), row.names = c(
    NA,
    -2L
  ), class = c("tbl_df", "tbl", "data.frame")), structure(list(
    .time = c(20, 100), .pred_survival = c(
      0.97502549742525,
      0.509525293895381
    ), .pred_lower = c(0.972295395585556, 0.498684513554291), .pred_upper = c(0.977535050974671, 0.52014429426975)
  ), row.names = c(
    NA,
    -2L
  ), class = c("tbl_df", "tbl", "data.frame")))), class = c(
    "tbl_df",
    "tbl", "data.frame"
  ), row.names = c(NA, -2L))

  expect_equal(pred, expected)
})

test_that("hazard prediction works", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  new_data <- data.frame(x = c("0", "1"))
  pred <- predict(mod, new_data = new_data, type = "hazard", eval_time = c(20, 100), interval = "credible")

  expected <- structure(list(.pred = list(structure(list(
    .time = c(20, 100),
    .pred_hazard = c(0.014685958026776, 0.0557037183479297),
    .pred_lower = c(0.0140749390839117, 0.0482821138874736),
    .pred_upper = c(0.0152389992297559, 0.0657579719365234)
  ), row.names = c(
    NA,
    -2L
  ), class = c("tbl_df", "tbl", "data.frame")), structure(list(
    .time = c(20, 100), .pred_hazard = c(
      0.00470160254288338,
      0.0121971535806942
    ), .pred_lower = c(
      0.00431229071012346,
      0.0113636203373507
    ), .pred_upper = c(0.005091475834385, 0.0130761077120619)
  ), row.names = c(NA, -2L), class = c("tbl_df", "tbl", "data.frame")))), class = c("tbl_df", "tbl", "data.frame"), row.names = c(
    NA,
    -2L
  ))

  expect_equal(pred, expected)
})

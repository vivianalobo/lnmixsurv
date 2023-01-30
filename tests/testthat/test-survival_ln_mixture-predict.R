test_that("survival prediction works", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  new_data <- data.frame(x = c("0", "1"))
  pred <- predict(mod, new_data, type = "survival", time = c(20, 100))

  expected <- structure(list(.pred = list(structure(list(
    .time = c(20, 100),
    .pred_survival = c(0.890693737049813, 0.00982157645355162)
  ), class = c("tbl_df", "tbl", "data.frame"), row.names = c(
    NA,
    -2L
  )), structure(list(.time = c(20, 100), .pred_survival = c(
    0.97502549742525,
    0.509525293895381
  )), class = c("tbl_df", "tbl", "data.frame"), row.names = c(
    NA,
    -2L
  )))), class = c("tbl_df", "tbl", "data.frame"), row.names = c(
    NA,
    -2L
  ))

  expect_equal(pred, expected)
})

test_that("hazard prediction works", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  new_data <- data.frame(x = c("0", "1"))
  pred <- predict(mod, new_data, type = "hazard", time = c(20, 100))

  expected <- structure(list(.pred = list(structure(list(
    .time = c(20, 100),
    .pred_hazard = c(0.014685958026776, 0.0557037183479297)
  ), class = c(
    "tbl_df",
    "tbl", "data.frame"
  ), row.names = c(NA, -2L)), structure(list(
    .time = c(20, 100), .pred_hazard = c(
      0.00470160254288338,
      0.0121971535806942
    )
  ), class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -2L)))), class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -2L))

  expect_equal(pred, expected)
})

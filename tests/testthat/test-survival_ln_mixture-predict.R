test_that("survival prediction works", {
  mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
  new_data <- data.frame(x = c("0", "1"))
  pred <- predict(mod, new_data, type = "survival", eval_time = c(20, 100), interval = "credible")

  expected <- structure(
    list(.pred = list(structure(list(.eval_time = c(
      20,
      100
    ), .pred_survival = c(0.890695166273467, 0.00959208655672782), .pred_lower = c(0.884838116393057, 0.00730165088523564), .pred_upper = c(
      0.896598543618103,
      0.0120560317734822
    )), row.names = c(NA, -2L), class = c(
      "tbl_df",
      "tbl", "data.frame"
    )), structure(list(
      .eval_time = c(20, 100),
      .pred_survival = c(0.975038970820135, 0.510000163508516),
      .pred_lower = c(0.972255921061828, 0.498098556543236), .pred_upper = c(
        0.977622037662899,
        0.521275512802539
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
  ), .pred_hazard = c(0.0146964470697921, 0.0568396736398598), .pred_lower = c(0.0140778684863756, 0.0496863207587222), .pred_upper = c(
    0.0152862604886661,
    0.0684146856583835
  )), row.names = c(NA, -2L), class = c(
    "tbl_df",
    "tbl", "data.frame"
  )), structure(list(
    .eval_time = c(20, 100),
    .pred_hazard = c(0.00470521245746177, 0.0122878523767787),
    .pred_lower = c(0.00434427274450936, 0.0114296201171799),
    .pred_upper = c(0.00509330799796965, 0.0132225119770675)
  ), row.names = c(
    NA,
    -2L
  ), class = c("tbl_df", "tbl", "data.frame")))), class = c(
    "tbl_df",
    "tbl", "data.frame"
  ), row.names = c(NA, -2L))

  expect_equal(pred, expected)
})

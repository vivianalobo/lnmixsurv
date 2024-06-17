# test_that("survival prediction works", {
#   mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
#   new_data <- data.frame(x = c("0", "1"))
#   pred <- predict(mod, new_data, type = "survival", eval_time = c(20, 100), interval = "credible")
#
#   expected <- structure(
#     list(.pred = list(structure(list(
#       .eval_time = c(
#         20,
#         100
#       ), .pred_survival = c(0.891151125, 0.009867837),
#       .pred_lower = c(0.885160087, 0.008194163),
#       .pred_upper = c(
#         0.89670372,
#         0.01204335
#       )
#     ), row.names = c(NA, -2L), class = c(
#       "tbl_df",
#       "tbl", "data.frame"
#     )), structure(list(
#       .eval_time = c(20, 100),
#       .pred_survival = c(0.9751848, 0.5098393),
#       .pred_lower = c(0.9721853, 0.4977577), .pred_upper = c(
#         0.9777597,
#         0.5211575
#       )
#     ), row.names = c(NA, -2L), class = c(
#       "tbl_df",
#       "tbl", "data.frame"
#     )))),
#     class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -2L)
#   )
#
#   expect_equal(pred, expected, tolerance = 10^-1)
# })
#
# test_that("hazard prediction works", {
#   mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
#   new_data <- data.frame(x = c("0", "1"))
#   pred <- predict(mod, new_data = new_data, type = "hazard", eval_time = c(20, 100), interval = "credible")
#
#   expected <- structure(list(.pred = list(structure(list(
#     .eval_time = c(
#       20,
#       100
#     ), .pred_hazard = c(0.01466505, 0.05542897),
#     .pred_lower = c(0.01404703, 0.04918021),
#     .pred_upper = c(
#       0.01525981,
#       0.06235824
#     )
#   ), row.names = c(NA, -2L), class = c(
#     "tbl_df",
#     "tbl", "data.frame"
#   )), structure(list(
#     .eval_time = c(20, 100),
#     .pred_hazard = c(0.004687204, 0.012180635),
#     .pred_lower = c(0.004296136, 0.011277945),
#     .pred_upper = c(0.005080489, 0.013036060)
#   ), row.names = c(
#     NA,
#     -2L
#   ), class = c("tbl_df", "tbl", "data.frame")))), class = c(
#     "tbl_df",
#     "tbl", "data.frame"
#   ), row.names = c(NA, -2L))
#
#   expect_equal(pred, expected, tolerance = 10^-1)
# })

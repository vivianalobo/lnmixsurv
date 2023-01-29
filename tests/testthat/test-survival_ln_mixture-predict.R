test_that("survival prediction works", {
    mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
    new_data <- data.frame(x = c("0", "1"))
    pred <- predict(mod, new_data, type = "survival", time = c(20, 100))

    expected <- structure(list(.pred = list(structure(list(
        .time = c(20, 100),
        .pred_survival = c(0.97787450088009, 0.0785202948828178)
    ), class = c(
        "tbl_df",
        "tbl", "data.frame"
    ), row.names = c(NA, -2L)), structure(list(
        .time = c(20, 100), .pred_survival = c(
            0.999769440516025,
            0.331594975084514
        )
    ),
    class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -2L)
    ))), class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -2L))

    expect_equal(pred, expected)
})

test_that("hazard prediction works", {
    mod <- readRDS(test_path("fixtures", "ln_fit_with_covariates.rds"))
    new_data <- data.frame(x = c("0", "1"))
    pred <- predict(mod, new_data, type = "hazard", time = c(20, 100))

    expected <- structure(list(.pred = list(structure(list(
        .time = c(20, 100),
        .pred_hazard = c(0.0077241222678161, 0.0272029126776698)
    ), class = c(
        "tbl_df",
        "tbl", "data.frame"
    ), row.names = c(NA, -2L)), structure(list(
        .time = c(20, 100), .pred_hazard = c(
            7.39023502584329e-05,
            0.00965847257456305
        )
    ),
    class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -2L)
    ))), class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -2L))

    expect_equal(pred, expected)
})

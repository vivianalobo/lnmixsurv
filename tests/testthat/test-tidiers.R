test_that("tiders work as expected", {
    mod <- new_survival_ln_fit_with_fixed_rng()
    obtained <- tidy(mod, conf.int = TRUE)
    expected <- structure(list(
        term = c(
            "beta_a[1]", "beta_a[2]", "beta_b[1]",
            "beta_b[2]"
        ), estimate = c(
            3.42640069809889, 0.488854437096543,
            4.04480760255212, 0.809520356761791
        ), std.error = c(
            0.0223196277628254,
            0.0191176391183816, 0.00635669062167443, 0.0100016430850216
        ),
        conf.low = c(
            3.40152555804949, 0.461572125886593, 4.03645285526598,
            0.79709270483321
        ), conf.high = c(
            3.45530477038196, 0.514435660330035,
            4.05284276064431, 0.822129664438226
        )
    ), row.names = c(
        NA,
        -4L
    ), class = c("draws_summary", "tbl_df", "tbl", "data.frame"))
    expect_equal(obtained, expected)
})

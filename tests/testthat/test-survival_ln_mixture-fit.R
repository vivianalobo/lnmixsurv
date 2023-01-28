
test_that("outcome must be Surv object", {
    expect_error(
        survival_ln_mixture(y ~ x, sim_data$data)
    )
})

test_that("only right censored data supported", {
    expect_error(
        survival_ln_mixture(survival::Surv(y, delta, type = "left") ~ x, sim_data$data)
    )
})

test_that("when using ~NULL, intercept must be TRUE", {
    expect_error(
        survival_ln_mixture(survival::Surv(y, delta) ~ NULL, sim_data$data, intercept = FALSE)
    )
})

test_that("survival_ln_mixture doesnt work with xy specification", {
    expect_error(
        survival_ln_mixture(sim_data$data$y, sim_data$data$x)
    )
})

test_that("survival_ln_mixture works with intercept only fit",{
    mod <- new_survival_ln_fit_with_fixed_rng_intercept_only()
    expect_equal(tidy(mod)$estimate, c(3.8243392, 4.6485868))
})

test_that("fit works as expected with simulated data", {
    mod <- new_survival_ln_fit_with_fixed_rng()
    post_summary <- posterior::summarise_draws(mod$posterior, estimate = stats::median, std.error = stats::mad)
    colnames(post_summary)[1] <- "term"
    post_tidy <- tidy(mod, effects = c("fixed", "auxiliary"))
    expected_result <- structure(list(term = c(
        "beta_a[1]", "beta_a[2]", "beta_b[1]",
        "beta_b[2]", "phi_a", "phi_b", "theta_a"
    ), estimate = c(
        3.42640069809889,
        0.488854437096543, 4.04480760255212, 0.809520356761791, 26.5081614971505,
        3.17737806452232, 0.505594778717908
    ), std.error = c(
        0.0223196277628254,
        0.0191176391183816, 0.00635669062167443, 0.0100016430850216,
        1.37813558584402, 0.120429036758056, 0.013495809879222
    )), row.names = c(
        NA,
        -7L
    ), class = c("draws_summary", "tbl_df", "tbl", "data.frame"))

    expect_equal(mod$nobs, 10000)
    expect_equal(post_summary, post_tidy)
    expect_equal(post_summary, expected_result)
})

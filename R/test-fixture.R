globalVariables("sim_data")

new_survival_ln_fit_with_fixed_rng <- function() {
    withr::with_rng_version(
        "4.2.2",
        withr::with_seed(
            1,
            survival_ln_mixture(survival::Surv(y, delta) ~ x, sim_data$data) # nolint: object_usage_linter.
        )
    )
}

new_survival_ln_fit_with_fixed_rng_intercept_only <- function() {
    withr::with_rng_version(
        "4.2.2",
        withr::with_seed(
            1,
            survival_ln_mixture(
                survival::Surv(y, delta) ~ NULL, sim_data$data, # nolint: object_usage_linter.
                intercept = TRUE
            )
        )
    )
}

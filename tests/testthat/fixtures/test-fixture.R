globalVariables("sim_data")

path <- "tests/testthat/fixtures/"

ln_fit_with_covariates <- withr::with_rng_version(
    "4.2.2",
    withr::with_seed(
        1,
        survival_ln_mixture(survival::Surv(y, delta) ~ x, sim_data$data) # nolint: object_usage_linter.
    )
)
saveRDS(ln_fit_with_covariates, paste0(path, "ln_fit_with_covariates.rds"))

ln_fit_with_intercept_only <- withr::with_rng_version(
    "4.2.2",
    withr::with_seed(
        1,
        survival_ln_mixture(
            survival::Surv(y, delta) ~ NULL, sim_data$data, # nolint: object_usage_linter.
            intercept = TRUE
        )
    )
)
saveRDS(ln_fit_with_intercept_only, paste0(path, "ln_fit_with_intercept_only.rds"))

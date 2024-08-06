globalVariables("sim_data")

path <- "tests/testthat/fixtures/"

ln_fit_with_covariates <- survival_ln_mixture(survival::Surv(y, delta) ~ x, sim_data$data, 
                                              starting_seed = 10, em_iter = 50)

saveRDS(ln_fit_with_covariates, paste0(path, "ln_fit_with_covariates.rds"))

ln_fit_with_intercept_only <- survival_ln_mixture(
  survival::Surv(y, delta) ~ NULL, sim_data$data, # nolint: object_usage_linter.
  intercept = TRUE, starting_seed = 10, em_iter = 50
)

saveRDS(ln_fit_with_intercept_only, paste0(path, "ln_fit_with_intercept_only.rds"))

new_survival_ln_mixture <- function(posterior, nobs, blueprint) {
    hardhat::new_model(posterior = posterior, nobs = nobs, blueprint = blueprint, class = "survival_ln_mixture")
}

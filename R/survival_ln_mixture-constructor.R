new_survival_ln_mixture <- function(posterior, blueprint) {
    hardhat::new_model(posterior = posterior, blueprint = blueprint, class = "survival_ln_mixture")
}

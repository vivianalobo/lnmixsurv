new_survival_ln_mixture <- function(posterior, nobs, predictors_name, blueprint) {
  hardhat::new_model(
    posterior = posterior,
    nobs = nobs,
    predictors_name = predictors_name,
    blueprint = blueprint,
    class = "survival_ln_mixture"
  )
}

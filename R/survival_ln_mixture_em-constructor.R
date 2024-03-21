new_survival_ln_mixture_em <- function(em_iterations,
                                       nobs,
                                       predictors_name,
                                       mixture_groups,
                                       blueprint) {
  hardhat::new_model(
    em_iterations = em_iterations,
    nobs = nobs,
    predictors_name = predictors_name,
    mixture_groups = mixture_groups,
    blueprint = blueprint,
    class = "survival_ln_mixture_em"
  )
}
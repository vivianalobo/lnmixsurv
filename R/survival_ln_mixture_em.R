#' Lognormal mixture model - EM algorithm for censored data
#'
#' `survival_ln_mixture_em()` fits a lognormal mixture model using the EM algorithm for censored data, as described in LOBO, Viviana GR; FONSECA, Thaís CO; ALVES, Mariane B. Lapse risk modeling in insurance: a Bayesian mixture approach. Annals of Actuarial Science, v. 18, n. 1, p. 126-151, 2024.
#'
#' @param y Log of the events ocurrence times.
#' 
#' @param delta Vector containing the observations status (0 for censored, 1 for event ocurrence)
#' 
#' @param X Design matrix for regression. The columns are the covariates (intercept + regression covariates) and the rows are the observations.
#'
#' @param iter A positive integer specifying the number of iterations for the EM algorithm.
#'
#' @param mixture_components number of mixture componentes, >= 2.
#'
#' @note Categorical predictors must be converted to factors before the fit,
#' otherwise the predictions will fail.
#'
#' @export
survival_ln_mixture_em <- function(y, delta, X, iter = 50, mixture_components = 2) {
  matrix_em_iter <- lognormal_mixture_em(iter, mixture_components, y, delta, X)
  
  n_regressors <- ncol(X)
  new_names <- NULL

  for(g in 1:mixture_components) {
    for(j in 1:3) {
      if(j == 1) {
        new_names <- c(new_names,
                       paste0('eta_', g))
      } else if ( j == 2) {
        for(k in 0:(n_regressors - 1)) {
          new_names <- c(new_names,
                         paste0('beta', k, '_', g))
        }
      } else {
        new_names <- c(new_names,
                       paste0('phi_', g))
      }
    }
  }

  colnames(matrix_em_iter) <- new_names
  
  matrix_em_iter <- matrix_em_iter |> 
    tibble::as_tibble() |> 
    dplyr::mutate(iter = 1:Niter)
  
  matrix_em_iter <- matrix_em_iter |> 
    tidyr::pivot_longer(1:(ncol(matrix_em_iter) - 1))
  
  return(matrix_em_iter)
}


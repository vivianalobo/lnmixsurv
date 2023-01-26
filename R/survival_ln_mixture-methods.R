#' Extract the number of observations from `survival_ln_mixture` fit.
#' 
#' Extract the number of observations used in a `survival_ln_mixture` fit.
#' 
#' @param model A fitted `survival_ln_mixture` object.
#' 
#' @return A single integer.
#' 
#' @export
nobs.survival_ln_mixture <- function(model, ...) {
  return(model$nobs)
}

extract_posterior <- function(model) {
  return(model$posterior)
}

extract_formula <- function(model) {
  # Trocar NULL por 1 para caso onde so tem intercepto
  formula <- gsub("NULL", "1", deparse(a$blueprint$formula))
  return(formula)
}

npredictors <- function(model) {
  return(ncol(model$blueprint$ptypes$predictors) + model$blueprint$intercept)
}

coef_names <- function(model) {
  intercepto <- ifelse(model$blueprint$intercept, "(Intercept)", "")
  cov <- colnames(model$blueprint$ptypes$predictors)
  return(c(intercepto, cov))
}

# Eu quero um método para extrair os coeficientes?
# acho que teria q sumarizar a posteriori junto com o ajuste do modelo
# não sei se é o ideal. Talvez seja melhor usar o tidy para extrair essas
# coisas mesmo.
# #' @export
# coef.survival_ln_mixture <- function(object, ...) {
#     post <- extract_posterior(object)
#     post <- posterior::merge_chains(post)
#     post <- posterior::subset_draws(post, variable = c("beta_a", "beta_b"))
#     return(as.vector(posterior::summarise_draws(post, median)$median))
# }

# se.survival_ln_mixture <- function(object, ...){
#     post <- extract_posterior(object)
#     post <- posterior::merge_chains(post)
#     post <- posterior::subset_draws(post, variable = c("beta_a", "beta_b"))
#     return(as.vector(posterior::summarise_draws(post, mad)$mad))
# }

# sigma.survival_ln_mixture <- function(object, ...){
#     post <- extract_posterior(object)
#     post <- posterior::merge_chains(post)
#     post <- posterior::subset_draws(post, variable = c("beta_a", "beta_b"))
#     return(as.vector(posterior::summarise_draws(post, mad)$mad))
# }

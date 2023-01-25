nobs <- function(model) {
    return(model$nobs)
}

extract_posterior <- function(model) {
    return(model$posterior)
}

extract_formula <- function(model) {
    return(deparse(model$blueprint$formula))
}

npredictors <- function(model) {
    return(ncol(model$blueprint$ptypes$predictors))
}

coef_names <- function(model) {
    return(colnames(model$blueprint$ptypes$predictors))
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
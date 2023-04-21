#' Lognormal mixture model
#'
#' `survival_ln_mixture()` fits a lognormal mixture model, as described in (referencia artigo viviana).
#' Colocar mais detalhes sobre o modelo (como a equacao) em details.
#'
#' @param formula A formula specifying the outcome terms on the left-hand side,
#' and the predictor terms on the right-hand side. The outcome must be a [survival::Surv]
#' object.
#'
#' @param data A __data frame__ containing both the predictors and the outcome.
#'
#' @param intercept A logical. Should an intercept be included in the processed data?
#'
#' @param iter A positive integer specifying the number of iterations for each chain (including warmup).
#'
#' @param warmup A positive integer specifying the number of warmup (aka burnin) iterations per chain.
#' The number of warmup iterations should be smaller than iter.
#'
#' @param thin A positive integer specifying the period for saving samples.
#'
#' @param chains A positive integer specifying the number of Markov chains.
#'
#' @param cores Ignored. Parallel runs are disabled.
#' 
#' @param numero_componentes number of mixture componentes. Currently, only accepts 2 or 3.
#'
#' @param ... Not currently used, but required for extensibility.
#'
#' @note Categorical predictos must be converted to factores before the fit,
#' otherwise the predictions will fail.
#'
#' @return
#'
#' A `survival_ln_mixture` object, which is a list with the following componentes:
#'
#' \item{posterior}{A [posterior::draws_matrix] with the posterior of the parameters of the model.}
#' \item{nobs}{A integer holding the number of observations used to generate the fit.}
#' \item{blueprint}{The blueprint component of the output of [hardhat::mold]}
#'
#'
#' @examples
#'
#' # Formula interface
#' library(survival)
#' set.seed(1)
#' mod <- survival_ln_mixture(Surv(time, status == 2) ~ NULL, lung, intercept = TRUE)
#'
#' @export
survival_ln_mixture <- function(formula, data, intercept = TRUE, iter = 1000, warmup = floor(iter / 10),
                                thin = 1, chains = 1, cores = 1, numero_componentes = 2, ...) {
  if (!numero_componentes %in% c(2,3)) stop("Only 2 or 3 componentes supported.") 
  rlang::check_dots_empty(...)
  UseMethod("survival_ln_mixture")
}

#' @export
#' @rdname survival_ln_mixture
survival_ln_mixture.default <- function(formula, ...) {
  stop("`survival_ln_mixture()` is not defined for a '", class(formula)[1], "'.", call. = FALSE)
}

# Formula method

#' @export
#' @rdname survival_ln_mixture
survival_ln_mixture.formula <- function(formula, data, intercept = TRUE, ...) {
  blueprint <- hardhat::default_formula_blueprint(intercept = intercept)
  processed <- hardhat::mold(formula, data, blueprint = blueprint)
  survival_ln_mixture_bridge(processed, ...)
}

# ------------------------------------------------------------------------------
# Bridge

survival_ln_mixture_bridge <- function(processed, ...) {
  predictors <- as.matrix(processed$predictors)
  outcome <- processed$outcome[[1]]

  if (!survival::is.Surv(outcome)) {
    rlang::abort("Response must be a survival object (created with survival::Surv)")
  }
  if (attr(outcome, "type") != "right") rlang::abort("Only right-censored data allowed")

  outcome_times <- outcome[, 1]
  outcome_status <- outcome[, 2]

  fit <- survival_ln_mixture_impl(predictors, outcome_times, outcome_status, ...)

  new_survival_ln_mixture(
    posterior = fit$posterior,
    nobs = fit$nobs,
    predictors_name = fit$predictors_name,
    blueprint = processed$blueprint
  )
}


# ------------------------------------------------------------------------------
# Implementation

survival_ln_mixture_impl <- function(predictors, outcome_times, outcome_status,
                                     iter = 1000, warmup = floor(iter / 10), thin = 1,
                                     chains = 1, cores = 1, numero_componentes = 2) {
  number_of_predictors <- ncol(predictors)
  if (number_of_predictors < 1) {
    rlang::abort(
      c(
        "The model must contain at least one predictor.",
        i = "When using outcome ~ NULL, intercept must be explicitly set to TRUE."
      )
    )
  }
  if (cores != 1) warning("Argumento cores ignorado, rodando cadeias sequencialmente.")
  
  fun <- switch(numero_componentes - 1,
    sequential_lognormal_mixture_gibbs_2_componentes,
    sequential_lognormal_mixture_gibbs_3_componentes
  )
  
  posterior_dist <- fun(predictors, outcome_times, outcome_status, iter, chains, 0)
  

  posterior_dist <- abind::abind(posterior_dist)

  grupos <- letters[1:numero_componentes]
  preds <- colnames(predictors)
  names_beta <- glue::glue_data(expand.grid(preds, grupos), "{Var1}_{Var2}")
  names_phi <- glue::glue("phi_{grupos}")
  names_theta <- glue::glue("theta_{grupos}")

  dimnames(posterior_dist)[[3]] <- c(names_beta, names_phi, names_theta)

  # Ajustar labels.
  theta <- apply(posterior_dist[, , names_theta, drop = FALSE], 3, stats::median)
  
  ordem <- order(theta, decreasing = TRUE)
  label_old <- dimnames(posterior_dist)[[3]]
  label_new <- c(
    glue::glue_data(expand.grid(preds, grupos[ordem]), "{Var1}_{Var2}"),
    names_phi[ordem],
    names_theta[ordem]
  )
  posterior_dist[,,label_old] <- posterior_dist[,,label_new]
  remover_menor_theta = -which(dimnames(posterior_dist)[[3]] == names_theta[numero_componentes])
  posterior_dist <- posterior_dist[,,remover_menor_theta]
  

  posterior_dist <- posterior::as_draws_matrix(posterior_dist)
  posterior_dist <- posterior::subset_draws(posterior_dist, iteration = seq(from = warmup + 1, to = iter))
  posterior_dist <- posterior::thin_draws(posterior_dist, thin = thin)

  list(posterior = posterior_dist, nobs = length(outcome_times), predictors_name = preds)
}

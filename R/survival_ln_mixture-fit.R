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
#' @param cores Number of cores to use when executing the chains in parallel. The number of
#' cores used is limited by the number of avaliable cores and the number of chains. The number
#' of available cores can be checked with `parallel::detectCores(logical = FALSE)`.
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
                                thin = 1, chains = 1, cores = 1, ...) {
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
    blueprint = processed$blueprint
  )
}


# ------------------------------------------------------------------------------
# Implementation

survival_ln_mixture_impl <- function(predictors, outcome_times, outcome_status,
                                     iter = 1000, warmup = floor(iter / 10), thin = 1,
                                     chains = 1, cores = 1) {
  if (chains > 1) rlang::warn("Check the posterior draws for label switch problem.")
  if (cores > 1) {
    posterior_dist <- parallel_lognormal_mixture_gibbs(
      predictors, outcome_times, outcome_status, iter, chains, cores, 0
    )
  } else {
    posterior_dist <- sequential_lognormal_mixture_gibbs(
      predictors, outcome_times, outcome_status, iter, chains, 0
    )
  }

  posterior_dist <- abind::abind(posterior_dist)

  grupos <- c("a", "b")
  preds <- seq_len(ncol(predictors))
  names_beta <- glue::glue_data(expand.grid(preds, grupos), "beta_{Var2}[{Var1}]")
  names_phi <- glue::glue("phi_{grupos}")

  dimnames(posterior_dist)[[length(dimnames(posterior_dist))]] <- c(names_beta, names_phi, "theta_a")

  posterior_dist <- posterior::as_draws_matrix(posterior_dist)
  posterior_dist <- posterior::subset_draws(posterior_dist, iteration = seq(from = warmup + 1, to = iter))
  posterior_dist <- posterior::thin_draws(posterior_dist, thin = thin)

  list(posterior = posterior_dist, nobs = length(outcome_times))
}

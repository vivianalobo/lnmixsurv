#' Fit a `survival_ln_mixture`
#'
#' `survival_ln_mixture()` fits a model.
#'
#' @param data When a __recipe__ or __formula__ is used, `data` is specified as:
#'
#'   * A __data frame__ containing both the predictors and the outcome.
#'   * The outcome must be a __survival::Surv__ objetc.
#'
#' @param formula A formula specifying the outcome terms on the left-hand side,
#' and the predictor terms on the right-hand side.
#'
#'   * The outcome must be a __survival::Surv__ object.
#'
#'
#' @param ... Not currently used, but required for extensibility.
#'
#' @return
#'
#' A `survival_ln_mixture` object.
#'
#' @examples
#'
#' # Formula interface
#' library(survival)
#' mod <- survival_ln_mixture(Surv(time, status == 2), lung)
#'
#'
#' # Recipes interface
#' library(recipes)
#' dados <- lung %>% mutate(outcome = Surv(time, status == 2))
#' rec <- recipe(outcome ~ sex, dados)
#' mod3 <- survival_ln_mixture(rec, dados)
#'
#' @export
survival_ln_mixture <- function(x, ...) {
    UseMethod("survival_ln_mixture")
}

#' @export
#' @rdname survival_ln_mixture
survival_ln_mixture.default <- function(x, ...) {
    stop("`survival_ln_mixture()` is not defined for a '", class(x)[1], "'.", call. = FALSE)
}

# Formula method

#' @export
#' @rdname survival_ln_mixture
survival_ln_mixture.formula <- function(formula, data, intercept = TRUE, ...) {
    blueprint <- hardhat::default_formula_blueprint(intercept = intercept)
    processed <- hardhat::mold(formula, data, blueprint = blueprint)
    survival_ln_mixture_bridge(processed, ...)
}

# Recipe method

#' @export
#' @rdname survival_ln_mixture
survival_ln_mixture.recipe <- function(x, data, intercept = TRUE, ...) {
    blueprint <- hardhat::default_recipe_blueprint(intercept = intercept)
    processed <- hardhat::mold(x, data, blueprint = blueprint)
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
        posterior_dist <- parallel_lognormal_mixture_gibbs_cpp(
            predictors, outcome_times, outcome_status, iter, chains, cores, 0
        )
    } else {
        posterior_dist <- sequential_lognormal_mixture_gibbs_cpp(
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

    list(posterior = posterior_dist)
}

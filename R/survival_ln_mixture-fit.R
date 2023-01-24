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
survival_ln_mixture.formula <- function(formula, data, ...) {
    processed <- hardhat::mold(formula, data)
    survival_ln_mixture_bridge(processed, ...)
}

# Recipe method

#' @export
#' @rdname survival_ln_mixture
survival_ln_mixture.recipe <- function(x, data, ...) {
    processed <- hardhat::mold(x, data)
    survival_ln_mixture_bridge(processed, ...)
}

# ------------------------------------------------------------------------------
# Bridge

survival_ln_mixture_bridge <- function(processed, ...) {
    predictors <- as.matrix(processed$predictors)
    outcome <- processed$outcome[[1]]

    if (!survival::is.Surv(outcome)) rlang::abort("Response must be a survival object (created with survival::Surv)")
    if (attr(outcome, "type") != "right") rlang::abort("Only right-censored data allowed")

    outcome_times <- outcome[, 1]
    outcome_status <- outcome[, 2]

    fit <- survival_ln_mixture_impl(predictors, outcome_times, outcome_status)

    new_survival_ln_mixture(
        posterior = fit$posterior,
        blueprint = processed$blueprint
    )
}


# ------------------------------------------------------------------------------
# Implementation

survival_ln_mixture_impl <- function(predictors, outcome_times, outcome_status) {
    posterior_dist <- lognormal_mixture_gibbs_cpp(outcome_times, predictors, outcome_status, 1000, 0)
    posterior_dist <- abind::abind(posterior_dist)

    grupos <- c("a", "b")
    preds <- seq_len(ncol(predictors))
    names_beta <- glue::glue_data(expand.grid(preds, grupos), "beta_{Var2}[{Var1}]")
    names_phi <- glue::glue("phi_{grupos}")

    dimnames(posterior_dist)[[length(dimnames(posterior_dist))]] <- c(names_beta, names_phi, "theta_a")

    posterior_dist <- posterior::as_draws_matrix(posterior_dist)

    list(posterior = posterior_dist)
}
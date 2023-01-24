#' Predict from a `survival_ln_mixture`
#'
#' @param object A `survival_ln_mixture` object.
#'
#' @param new_data A data frame or matrix of new predictors.
#'
#' @param type A single character. The type of predictions to generate.
#' Valid options are:
#'
#' - `"time"` for the survival time.
#' - `"survival"` for the survival probability.
#' - `"hazard"` for the hazard.
#'
#' @param ... Not used, but required for extensibility.
#'
#' @return
#'
#' A tibble of predictions. The number of rows in the tibble is guaranteed
#' to be the same as the number of rows in `new_data`.
#'
#' @examples
#' train <- mtcars[1:20, ]
#' test <- mtcars[21:32, -1]
#'
#' # Fit
#' mod <- survival_ln_mixture(mpg ~ cyl + log(drat), train)
#'
#' # Predict, with preprocessing
#' predict(mod, test)
#'
#' @export
predict.survival_ln_mixture <- function(object, new_data, type = "survival", ...) {
    forged <- hardhat::forge(new_data, object$blueprint)
    rlang::arg_match(type, valid_survival_ln_mixture_predict_types())
    predict_survival_ln_mixture_bridge(type, object, forged$predictors, ...)
}

valid_survival_ln_mixture_predict_types <- function() {
    c("time", "survival", "hazard")
}

# ------------------------------------------------------------------------------
# Bridge

predict_survival_ln_mixture_bridge <- function(type, model, predictors, ...) {
    predictors <- as.matrix(predictors)

    predict_function <- get_survival_ln_mixture_predict_function(type)
    predictions <- predict_function(model, predictors, ...)

    hardhat::validate_prediction_size(predictions, predictors)

    predictions
}

get_survival_ln_mixture_predict_function <- function(type) {
    switch(type,
        time = predict_survival_ln_mixture_time,
        survival = predict_survival_ln_mixture_survival,
        hazard = predict_survival_ln_mixture_hazard
    )
}

# ------------------------------------------------------------------------------
# Implementation

extract_surv_haz <- function(model, predictors, time, type = "survival") {
    rlang::arg_match(type, c("survival", "hazard"))

    fun <- switch(type,
        survival = sob_lognormal_mix,
        hazard = falha_lognormal_mix
    )

    post <- model$posterior

    qntd_chains <- psoterior::nchains(post)
    if (qntd_chains > 1) {
        rlang::abort("Merge the chains before predicting.")
    }

    qntd_iteracoes <- posterior::niterations(post)
    beta_a <- posterior::subset_draws(post, "beta_a")
    beta_b <- posterior::subset_draws(post, "beta_b")
    phi_a <- posterior::subset_draws(post, "phi_a")
    phi_b <- posterior::subset_draws(post, "phi_b")
    theta_a <- posterior::subset_draws(post, "theta_a")
    m_a <- beta_a %*% t(predictors)
    m_b <- beta_b %*% t(predictors)
    sigma_a <- sqrt(1 / phi_a)
    sigma_b <- sqrt(1 / phi_b)

    surv_haz <- list()
    for (i in seq_len(ncol(m_a))) {
        surv_haz[[i]] <- vapply(
            time, function(t) fun(t, m_a[, i], m_b[, i], sigma_a, sigma_b, theta_a), numeric(qntd_iteracoes)
        )
    }
    predictions <- lapply(surv_haz, function(x) apply(x, 2, stats::median))

    pred_name <- paste0(".pred_", type) # nolint: object_usage_linter.

    tibble::tibble(.pred = purrr::map(predictions, ~ tibble::tibble(.time = time, !!pred_name := .x)))
}

predict_survival_ln_mixture_time <- function(model, predictors) {
    rlang::abort("Not implemented")
    # predictions <- rep(1L, times = nrow(predictors))
    # hardhat::spruce_numeric(predictions)
}

predict_survival_ln_mixture_survival <- function(model, predictors, time) {
    extract_surv_haz(model, predictors, time, "survival")
}

predict_survival_ln_mixture_hazard <- function(model, predictors, time) {
    extract_surv_haz(model, predictors, time, "hazard")
}

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

    # hardhat::validate_prediction_size(predictions, predictors)

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

predict_survival_ln_mixture_time <- function(model, predictors) {
    predictions <- rep(1L, times = nrow(predictors))
    hardhat::spruce_numeric(predictions)
}

predict_survival_ln_mixture_survival <- function(model, predictors, time) {
    post <- model$posterior
    qntd_iteracoes <- nrow(post[[1]])
    sob <- matrix(NA, nrow = qntd_iteracoes, ncol = length(time))
    ma <- post[[1, 1]] %*% t(predictors)
    mb <- post[[2, 1]] %*% t(predictors)
    sigmaa <- sqrt(1 / post[[3, 1]])
    sigmab <- sqrt(1 / post[[4, 1]])
    theta <- post[[5, 1]]
    for (iteration in seq(qntd_iteracoes)) {
        sob[iteration, ] <- sob_lognormal_mix(
            time, ma[iteration], mb[iteration], sigmaa[iteration], sigmab[iteration], theta[iteration]
        )
    }
    predictions <- apply(sob, 2, stats::median)
    tibble::tibble(.time = time, .pred_survival = predictions)
    # hardhat::spruce_numeric(predictions)
}

predict_survival_ln_mixture_hazard <- function(model, predictors) {
    predictions <- rep(3L, times = nrow(predictors))
    hardhat::spruce_numeric(predictions)
}

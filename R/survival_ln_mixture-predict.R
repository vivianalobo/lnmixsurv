#' Predict from a Lognormal Mixture Model
#'
#' @param object A `survival_ln_mixture` object.
#'
#' @param new_data A data frame or matrix of new predictors.
#'
#' @param type A single character. The type of predictions to generate.
#' Valid options are:
#'
#' - `"time"` for the survival time. **not implmeented**
#' - `"survival"` for the survival probability.
#' - `"hazard"` for the hazard.
#'
#' @param eval_time For type = "hazard" or type = "survival", the times for the distribution.
#'
#' @param interval should interval estimates be added? Options are "none" and "credible".
#'
#' @param level the tail area of the intervals. Default value is 0.95.
#'
#' @param ... Not used, but required for extensibility.
#'
#' @note Categorical predictors must be converted to factors before the fit,
#' otherwise the predictions will fail.
#' @return
#'
#' A tibble of predictions. The number of rows in the tibble is guaranteed
#' to be the same as the number of rows in `new_data`.
#'
#' @examples
#'
#' # Categorical variables must be converted to factor before the fit.
#' require(survival)
#' set.seed(1)
#' mod <- survival_ln_mixture(Surv(time, status == 2) ~ factor(sex), lung, intercept = TRUE)
#' # Would result in error
#' \dontrun{
#' predict(mod, data.frame(sex = 1), eval_time = 100)
#' }
#'
#' # Correct way
#' lung$sex <- factor(lung$sex)
#' set.seed(1)
#' mod2 <- survival_ln_mixture(Surv(time, status == 2) ~ sex, lung, intercept = TRUE)
#' # Note: the categorical predictors must be character.
#' predict(mod2, data.frame(sex = "1"), eval_time = 100)
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

predict_survival_ln_mixture_time <- function(model, predictors) {
  rlang::abort("Not implemented")
}

predict_survival_ln_mixture_survival <- function(model, predictors, eval_time, interval = "none", level = 0.95) {
  extract_surv_haz(model, predictors, eval_time, interval, level, "survival")
}

predict_survival_ln_mixture_hazard <- function(model, predictors, eval_time, interval = "none", level = 0.95) {
  extract_surv_haz(model, predictors, eval_time, interval, level, "hazard")
}

extract_surv_haz <- function(model, predictors, eval_time, interval = "none", level = 0.95, type = "survival") {
  rlang::arg_match(type, c("survival", "hazard"))
  rlang::arg_match(interval, c("none", "credible"))

  fun <- switch(type,
    survival = sob_lognormal_mix,
    hazard = falha_lognormal_mix
  )

  post <- model$posterior

  qntd_chains <- posterior::nchains(post)
  if (qntd_chains > 1) {
    post <- posterior::merge_chains(post)
  }
  qntd_iteracoes <- posterior::niterations(post)
  beta <- lapply(model$mixture_groups, function(x) posterior::subset_draws(post, glue::glue("{model$predictors_name}_{x}")))
  phi <- posterior::subset_draws(post, "phi", regex = TRUE)
  theta <- posterior::subset_draws(post, "theta", regex = TRUE)
  m <- lapply(beta, function(x) x %*% t(predictors))
  sigma <- sqrt(1 / phi)

  surv_haz <- list()
  for (i in seq_len(length(model$predictors_name))) {
    surv_haz[[i]] <- vapply(
      eval_time, function(t) fun(t, lapply(m, function(x) x[,i]), sigma, theta), numeric(qntd_iteracoes)
    )
  }
  predictions <- lapply(surv_haz, function(x) apply(x, 2, stats::median))
  pred_name <- paste0(".pred_", type) # nolint: object_usage_linter.
  pred <- purrr::map(predictions, ~ tibble::tibble(.time = eval_time, !!pred_name := .x))

  if (interval == "credible") {
    lower <- lapply(surv_haz, function(x) apply(x, 2, stats::quantile, probs = 1 - level))
    upper <- lapply(surv_haz, function(x) apply(x, 2, stats::quantile, level))
    pred <- purrr::pmap(list(pred, lower, upper), function(x, y, z) {
      x$.pred_lower <- y
      x$.pred_upper <- z
      return(x)
    })
  }

  tibble::tibble(.pred = pred)
}

sob_lognormal_mix <- function(t, m, sigma, theta) {
  componentes <- vapply(seq_len(length(m) - 1), function(x) sob_lognormal(t, m[[x]], sigma[,x]) * theta[, x], numeric(nrow(sigma)))
  componentes <- cbind(componentes, sob_lognormal(t, m[[length(m)]], sigma[,length(m)]) * (1 - apply(theta, 1, sum)))
  s <- apply(componentes, 1, sum)
  return(s)
}

falha_lognormal_mix <- function(t, m, sigma, theta) {
  sob_mix <- sob_lognormal_mix(t, m, sigma, theta)
  componentes <- vapply(seq_len(length(m) - 1), function(x) stats::dlnorm(t, m[[x]], sigma[,x]) * theta[, x], numeric(nrow(sigma)))
  componentes <- cbind(componentes, stats::dlnorm(t, m[[length(m)]], sigma[,length(m)]) * (1 - apply(theta, 1, sum)))
  dlnorm_mix <- apply(componentes, 1, sum)
  return(dlnorm_mix / sob_mix)
}

sob_lognormal <- function(t, m, sigma) {
  s <- stats::pnorm((-log(t) + m) / sigma)
  return(s)
}

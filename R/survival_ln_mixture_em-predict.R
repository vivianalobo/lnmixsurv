#' Predict from a lognormal_em Mixture Model fitted using EM algorithm.
#'
#' @param object A `survival_ln_mixture_em` object.
#'
#' @param new_data A data frame or matrix of new predictors.
#'
#' @param type A single character. The type of predictions to generate.
#' Valid options are:
#'
#' - `"survival"` for the survival probability.
#' - `"hazard"` for the hazard theoretical hazard.
#'
#' @param eval_time For type = "hazard" or type = "survival", the times for the distribution.
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
#' @export
predict.survival_ln_mixture_em <- function(object, new_data, type,
                                           eval_time, ...) {
  forged <- hardhat::forge(new_data, object$blueprint)
  rlang::arg_match(type, valid_survival_ln_mixture_em_predict_types())

  predict_survival_ln_mixture_em_bridge(
    type, object, forged$predictors,
    eval_time, ...
  )
}

valid_survival_ln_mixture_em_predict_types <- function() {
  c("survival", "hazard")
}

# ------------------------------------------------------------------------------
# Bridge
predict_survival_ln_mixture_em_bridge <- function(type, model, predictors,
                                                  eval_time, ...) {
  predictors <- as.matrix(predictors)
  predict_function <- get_survival_ln_mixture_em_predict_function(type)
  predictions <- predict_function(model, predictors, eval_time, ...)
  predictions
}

get_survival_ln_mixture_em_predict_function <- function(type) {
  switch(type,
    survival = predict_survival_ln_mixture_em_survival,
    hazard = predict_survival_ln_mixture_em_hazard
  )
}

# ------------------------------------------------------------------------------
# Implementation
predict_survival_ln_mixture_em_time <- function(model, predictors, eval_time) {
  rlang::abort("Not implemented")
}

predict_survival_ln_mixture_em_survival <- function(model, predictors, eval_time) {
  extract_surv_haz_em(model, predictors, eval_time, "survival")
}

predict_survival_ln_mixture_em_hazard <- function(model, predictors, eval_time) {
  extract_surv_haz_em(model, predictors, eval_time, "hazard")
}

extract_surv_haz_em <- function(model, predictors, eval_time, type = "survival") {
  rlang::arg_match(type, c("survival", "hazard"))

  last_row <- model$em_iterations[nrow(model$em_iterations), -ncol(model$em_iterations)]

  beta <- matrix(
    as.numeric(last_row[
      !startsWith(names(last_row), "eta") & !(startsWith(names(last_row), "phi"))
    ]),
    ncol = model$mixture_groups
  )

  phi <- as.numeric(last_row[startsWith(names(last_row), "phi")])
  eta <- as.numeric(last_row[startsWith(names(last_row), "eta")])
  
  sigma <- 1 / sqrt(phi)

  m <- apply(beta,
    MARGIN = 2,
    FUN = function(x) predictors %*% x
  )

  if (type == "survival") {
    out <- list()
    if (nrow(predictors) > 1) {
      for (r in 1:nrow(predictors)) {
        out_r <- tibble::tibble(
          .eval_time = eval_time,
          .pred_survival = NA
        )
        
        out_r$.pred_survival <- as.numeric(predict_survival_em_cpp(eval_time, m, sigma, eta, r))

        out[[r]] <- out_r |>
          dplyr::bind_cols(multiply_rows(predictors[r, ], length(eval_time)))
      }
    } else {
      out_r <- tibble::tibble(
        .eval_time = eval_time,
        .pred_survival = NA
      )
      
      out_r$.pred_survival <- as.numeric(
        predict_survival_em_cpp(eval_time, t(as.matrix(m)), sigma, eta, 1))

      out[[1]] <- out_r |>
        dplyr::bind_cols(multiply_rows(predictors[1, ], length(eval_time)))
    }
  } else if (type == "hazard") {
    out <- list()
    if (nrow(predictors) > 1) {
      for (r in 1:nrow(predictors)) {
        out_r <- tibble::tibble(
          .eval_time = eval_time,
          .pred_hazard = NA
        )
        
        out_r$.pred_hazard <- as.numeric(predict_hazard_em_cpp(eval_time, m, sigma, eta, r))

        out[[r]] <- out_r |>
          dplyr::bind_cols(multiply_rows(predictors[r, ], length(eval_time)))
      }
    } else {
      out_r <- tibble::tibble(
        .eval_time = eval_time,
        .pred_hazard = NA
      )
      
      out_r$.pred_hazard <- as.numeric(predict_hazard_em_cpp(eval_time, t(as.matrix(m)), sigma, eta, 1))

      out[[1]] <- out_r |> 
        dplyr::bind_cols(multiply_rows(predictors[1, ], length(eval_time)))
    }
  }

  return(tibble::tibble(.pred = out))
}

multiply_rows <- function(x, n) {
  out <- tibble::as_tibble(as.data.frame(t(x)))

  if (n > 1) {
    for (i in 1:(n - 1)) {
      out <- out |>
        dplyr::slice(1) |>
        dplyr::bind_rows(out)
    }
  }

  return(out)
}

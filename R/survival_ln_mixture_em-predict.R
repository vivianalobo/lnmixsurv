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

  last_row <- model$em_iterations[nrow(model$em_iterations), ]

  beta <- matrix(as.numeric(last_row[startsWith(names(last_row), "beta")]),
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

        for (i in 1:length(eval_time)) {
          t <- eval_time[i]
          out_r$.pred_survival[i] <- sob_lognormal_em_mix(t, m[r, ], sigma, eta)
        }

        tib_predictors <- tibble::as_tibble(matrix(rep(predictors[r, ], times = length(eval_time)), ncol = 2, byrow = TRUE), .name_repair = "minimal")

        names(tib_predictors) <- colnames(predictors)
        out_r <- dplyr::bind_cols(
          out_r,
          tib_predictors
        )

        out[[r]] <- out_r
      }
    } else {
      out_r <- tibble::tibble(
        .eval_time = eval_time,
        .pred_survival = NA
      )

      for (i in 1:length(eval_time)) {
        t <- eval_time[i]
        out_r$.pred_survival[i] <- sob_lognormal_em_mix(t, m, sigma, eta)
      }

      tib_predictors <- tibble::as_tibble(matrix(rep(predictors[1, ], times = length(eval_time)), ncol = 2, byrow = TRUE), .name_repair = "minimal")

      names(tib_predictors) <- colnames(predictors)
      out_r <- dplyr::bind_cols(
        out_r,
        tib_predictors
      )

      out[[1]] <- out_r
    }
  } else if (type == "hazard") {
    out <- list()
    if (nrow(predictors) > 1) {
      for (r in 1:nrow(predictors)) {
        out_r <- tibble::tibble(
          .eval_time = eval_time,
          .pred_hazard = NA
        )

        for (i in 1:length(eval_time)) {
          t <- eval_time[i]
          out_r$.pred_hazard[i] <- falha_lognormal_em_mix(t, m[r, ], sigma, eta)
        }

        tib_predictors <- tibble::as_tibble(matrix(rep(predictors[r, ], times = length(eval_time)), ncol = 2, byrow = TRUE), .name_repair = "minimal")

        names(tib_predictors) <- colnames(predictors)
        out_r <- dplyr::bind_cols(
          out_r,
          tib_predictors
        )

        out[[r]] <- out_r
      }
    } else {
      out_r <- tibble::tibble(
        .eval_time = eval_time,
        .pred_hazard = NA
      )

      for (i in 1:length(eval_time)) {
        t <- eval_time[i]
        out_r$.pred_hazard[i] <- falha_lognormal_em_mix(t, m, sigma, eta)
      }

      tib_predictors <- tibble::as_tibble(matrix(rep(predictors[1, ], times = length(eval_time)), ncol = 2, byrow = TRUE), .name_repair = "minimal")

      names(tib_predictors) <- colnames(predictors)
      out_r <- dplyr::bind_cols(
        out_r,
        tib_predictors
      )

      out[[1]] <- out_r
    }
  }

  return(tibble::tibble(.pred = out))
}

sob_lognormal_em <- function(t, m, sigma) {
  stats::pnorm((-log(t) + m) / sigma)
}

sob_lognormal_em_mix <- function(t, m, sigma, eta) {
  out <- 0
  for (g in 1:length(m)) {
    out <- out + eta[g] * sob_lognormal_em(t, m[g], sigma[g])
  }
  return(out)
}

falha_lognormal_em_mix <- function(t, m, sigma, eta) {
  sob_mix <- sob_lognormal_em_mix(t, m, sigma, eta)
  dlnorm_mix <- 0
  for (g in 1:length(m)) {
    dlnorm_mix <- dlnorm_mix + eta[g] * stats::dlnorm(t, m[g], sigma[g])
  }
  return(dlnorm_mix / sob_mix)
}

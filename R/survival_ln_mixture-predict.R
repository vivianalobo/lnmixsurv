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
#' predict(mod, data.frame(sex = 1), type = "survival", eval_time = 100)
#' }
#'
#' # Correct way
#' lung$sex <- factor(lung$sex)
#' set.seed(1)
#' mod2 <- survival_ln_mixture(Surv(time, status == 2) ~ sex, lung, intercept = TRUE)
#' # Note: the categorical predictors must be character.
#' predict(mod2, data.frame(sex = "1"), type = "survival", eval_time = 100)
#'
#' @export
predict.survival_ln_mixture <- function(object, new_data, type, eval_time, interval = "none", level = 0.95, ...) {
  forged <- hardhat::forge(new_data, object$blueprint)
  rlang::arg_match(type, valid_survival_ln_mixture_predict_types())

  predict_survival_ln_mixture_bridge(type, object, forged$predictors, eval_time, interval, level, ...)
}

valid_survival_ln_mixture_predict_types <- function() {
  c("time", "survival", "hazard")
}

# ------------------------------------------------------------------------------
# Bridge

predict_survival_ln_mixture_bridge <- function(type, model, predictors, eval_time, interval, level, ...) {
  predictors <- as.matrix(predictors)

  predict_function <- get_survival_ln_mixture_predict_function(type)
  predictions <- predict_function(model, predictors, eval_time, interval, level, ...)

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

predict_survival_ln_mixture_time <- function(model, predictors, eval_time, interval, level) {
  rlang::abort("Not implemented")
}

predict_survival_ln_mixture_survival <- function(model, predictors, eval_time, interval, level) {
  extract_surv_haz(model, predictors, eval_time, interval, level, "survival")
}

predict_survival_ln_mixture_hazard <- function(model, predictors, eval_time, interval, level) {
  extract_surv_haz(model, predictors, eval_time, interval, level, "hazard")
}

extract_surv_haz <- function(model, predictors, eval_time, interval = "none",
                             level = 0.95, type = "survival") {
  rlang::arg_match(type, c("survival", "hazard"))
  rlang::arg_match(interval, c("none", "credible"))

  post <- model$posterior
  qntd_chains <- posterior::nchains(post)

  if (type == "survival") {
    out <- list()
    if (nrow(predictors) > 1) {
      for (r in 1:nrow(predictors)) {
        out_r <- tibble::tibble(
          .eval_time = rep(eval_time, qntd_chains)
        )
        preds <- NULL

        for (chain in 1:qntd_chains) {
          post_chain <- posterior::subset_draws(post, chain = chain)
          beta <- lapply(model$mixture_groups, function(x) {
            names <- paste0(model$predictors_name, "_", x)
            return(posterior::subset_draws(post_chain, names))
          })

          phi <- posterior::subset_draws(post_chain, "phi", regex = TRUE)
          eta <- posterior::subset_draws(post_chain, "eta", regex = TRUE)
          sigma <- sqrt(1 / phi)

          preds_chain <- as.data.frame(predict_survival_gibbs_cpp(
            eval_time, predictors[r, ],
            beta, sigma, eta,
            interval == "credible", level
          ))
          preds <- dplyr::bind_rows(preds, dplyr::bind_cols(preds_chain, tibble::tibble(chain = factor(chain))))
        }

        if (interval == "credible") {
          colnames(preds) <- c(".pred_survival", ".pred_lower", ".pred_upper", "chain")
        } else {
          colnames(preds) <- c(".pred_survival", "chain")
        }

        out_r <- dplyr::bind_cols(out_r, preds)

        out[[r]] <- out_r |>
          dplyr::bind_cols(multiply_rows(predictors[r, ], length(rep(eval_time, qntd_chains))))
      }
    } else {
      out_r <- tibble::tibble(
        .eval_time = rep(eval_time, qntd_chains)
      )

      preds <- NULL

      for (chain in 1:qntd_chains) {
        post_chain <- posterior::subset_draws(post, chain = chain)
        beta <- lapply(model$mixture_groups, function(x) {
          names <- paste0(model$predictors_name, "_", x)
          return(posterior::subset_draws(post_chain, names))
        })

        phi <- posterior::subset_draws(post_chain, "phi", regex = TRUE)
        eta <- posterior::subset_draws(post_chain, "eta", regex = TRUE)
        sigma <- sqrt(1 / phi)

        preds_chain <- as.data.frame(predict_survival_gibbs_cpp(
          eval_time, predictors[1, ],
          beta, sigma, eta,
          interval == "credible", level
        ))
        preds <- dplyr::bind_rows(preds, dplyr::bind_cols(preds_chain, tibble::tibble(chain = factor(chain))))
      }

      if (interval == "credible") {
        colnames(preds) <- c(".pred_survival", ".pred_lower", ".pred_upper", "chain")
      } else {
        colnames(preds) <- c(".pred_survival", "chain")
      }

      out_r <- dplyr::bind_cols(out_r, preds)

      out[[1]] <- out_r |>
        dplyr::bind_cols(multiply_rows(predictors[1, ], length(rep(eval_time, qntd_chains))))
    }
  } else {
    out <- list()
    if (nrow(predictors) > 1) {
      for (r in 1:nrow(predictors)) {
        out_r <- tibble::tibble(
          .eval_time = rep(eval_time, qntd_chains)
        )

        preds <- NULL

        for (chain in 1:qntd_chains) {
          post_chain <- posterior::subset_draws(post, chain = chain)
          beta <- lapply(model$mixture_groups, function(x) {
            names <- paste0(model$predictors_name, "_", x)
            return(posterior::subset_draws(post_chain, names))
          })

          phi <- posterior::subset_draws(post_chain, "phi", regex = TRUE)
          eta <- posterior::subset_draws(post_chain, "eta", regex = TRUE)
          sigma <- sqrt(1 / phi)

          preds_chain <- as.data.frame(predict_hazard_gibbs_cpp(
            eval_time, predictors[r, ],
            beta, sigma, eta,
            interval == "credible", level
          ))
          preds <- dplyr::bind_rows(preds, dplyr::bind_cols(preds_chain, tibble::tibble(chain = factor(chain))))
        }

        if (interval == "credible") {
          colnames(preds) <- c(".pred_hazard", ".pred_lower", ".pred_upper", "chain")
        } else {
          colnames(preds) <- c(".pred_hazard", "chain")
        }

        out_r <- dplyr::bind_cols(out_r, preds)

        out[[r]] <- out_r |>
          dplyr::bind_cols(multiply_rows(predictors[r, ], length(rep(eval_time, qntd_chains))))
      }
    } else {
      out_r <- tibble::tibble(
        .eval_time = rep(eval_time, qntd_chains)
      )

      preds <- NULL

      for (chain in 1:qntd_chains) {
        post_chain <- posterior::subset_draws(post, chain = chain)
        beta <- lapply(model$mixture_groups, function(x) {
          names <- paste0(model$predictors_name, "_", x)
          return(posterior::subset_draws(post_chain, names))
        })

        phi <- posterior::subset_draws(post_chain, "phi", regex = TRUE)
        eta <- posterior::subset_draws(post_chain, "eta", regex = TRUE)
        sigma <- sqrt(1 / phi)

        preds_chain <- as.data.frame(predict_survival_gibbs_cpp(
          eval_time, predictors[1, ],
          beta, sigma, eta,
          interval == "credible", level
        ))
        preds <- dplyr::bind_rows(preds, dplyr::bind_cols(preds_chain, tibble::tibble(chain = factor(chain))))
      }

      if (interval == "credible") {
        colnames(preds) <- c(".pred_hazard", ".pred_lower", ".pred_upper", "chain")
      } else {
        colnames(preds) <- c(".pred_hazard", "chain")
      }

      out_r <- dplyr::bind_cols(out_r, preds)

      out[[1]] <- out_r |>
        dplyr::bind_cols(multiply_rows(predictors[1, ], length(rep(eval_time, qntd_chains))))
    }
  }

  tibble_out <- tibble::tibble(.pred = out)

  form <- model$blueprint$formula
  form <- stats::as.formula(paste0(paste(form)[[2]], " ~ ", paste(form)[[3]]))
  vars <- strsplit(paste(form)[[3]], " \\+ ")[[1]]

  if (!model$blueprint$intercept) {
    vars <- vars[vars != "0"]
  }

  for (k in 1:length(tibble_out$.pred)) {
    preds <- tibble_out$.pred[[k]]

    if (all(vars != "NULL")) {
      data <- model$data

      dict <- tibble::tibble(
        key = NA,
        match = NA,
        category = NA
      )

      predictors_data <- data |>
        dplyr::select(dplyr::all_of(vars))

      if (length(vars) == 1) {
        predictors_check_levels <- predictors_data
        names(predictors_check_levels) <- "x"
        if (length(levels(droplevels(predictors_check_levels$x))) == 1) {
          return(preds)
        }
      }

      time <- NA
      estimate <- NA
      .eval_time <- NA
      .pred_survival <- NA
      .pred_hazard <- NA

      for (c in 1:ncol(predictors_data)) {
        if (is.factor(predictors_data |>
          dplyr::select(dplyr::all_of(c)) |>
          dplyr::pull())) {
          key <- names(predictors_data)[c]
          var <- predictors_data |>
            dplyr::select(dplyr::all_of(c)) |>
            dplyr::pull()

          for (i in levels(var)) {
            match <- paste0(key, i)
            category <- i
            dict <- dplyr::bind_rows(dict, tibble::tibble(
              key = key,
              match = match,
              category = category
            ))
          }
        }
      }

      dict <- dict |>
        dplyr::slice(-1)

      strata_preds <- rep(NA, nrow(preds))

      for (i in 1:nrow(preds)) {
        if (model$blueprint$intercept) {
          phrase <- NULL
          next_index <- 0

          for (c in (5 + 2 * as.numeric(interval == "credible")):ncol(preds)) {
            if (c >= next_index) {
              match <- names(preds)[c]
              key_c <- dict$key[dict$match == match]
              if (match %in% dict$match) {
                dict_key_c <- dict |>
                  dplyr::filter(key == key_c)

                preds_key_c <- preds |>
                  dplyr::select(dplyr::starts_with(key_c)) |>
                  dplyr::slice(i)

                index_1 <- which(preds_key_c == 1)
                if (identical(index_1, integer(0))) {
                  category <- dict_key_c$category[which(!(dict_key_c$match %in% names(preds_key_c)))]
                } else {
                  category <- dict_key_c$category[which(dict_key_c$match == names(preds_key_c)[index_1])]
                }

                if (is.null(phrase)) {
                  phrase <- paste0(phrase, key_c, "=", category)
                } else {
                  phrase <- paste0(phrase, ", ", key_c, "=", category)
                }

                for (m in c:ncol(preds)) {
                  if (m == ncol(preds) & names(preds)[m] %in% dict_key_c$match) {
                    next_index <- ncol(preds) + 1
                  } else {
                    if (names(preds)[m] %in% dict_key_c$match) {
                      next
                    } else {
                      next_index <- m
                      break
                    }
                  }
                }
              } else {
                if (is.null(phrase)) {
                  phrase <- paste0(phrase, match, "=", preds[i, c])
                } else {
                  phrase <- paste0(phrase, ", ", match, "=", preds[i, c])
                }
              }
            }
          }
          strata_preds[i] <- phrase
        } else {
          phrase <- NULL
          next_index <- 0
          if (length(unique(dict$key)) == 1) {
            cols_remove <- dict |>
              dplyr::group_by(key) |>
              dplyr::slice(-dplyr::n()) |>
              dplyr::pull(match)
          } else {
            cols_remove <- NULL
            for (ky in unique(dict$key)) {
              if (sum(startsWith(names(preds), ky)) ==
                nrow(dict |> dplyr::filter(key == ky))) {
                cols_remove <- c(cols_remove, dict$match[dict$key == ky][1])
              }
            }
          }

          preds_modified <- preds |>
            dplyr::select(-dplyr::all_of(cols_remove))

          for (c in (4 + 2 * as.logical(interval == "credible")):ncol(preds_modified)) {
            if (c >= next_index) {
              match <- names(preds_modified)[c]
              key_c <- dict$key[dict$match == match]
              if (match %in% dict$match) {
                dict_key_c <- dict |>
                  dplyr::filter(key == key_c)

                preds_key_c <- preds_modified |>
                  dplyr::select(dplyr::starts_with(key_c)) |>
                  dplyr::slice(i)

                index_1 <- which(preds_key_c == 1)
                if (identical(index_1, integer(0))) {
                  category <- dict_key_c$category[which(!(dict_key_c$match %in% names(preds_key_c)))]
                } else {
                  category <- dict_key_c$category[which(dict_key_c$match == names(preds_key_c)[index_1])]
                }

                if (is.null(phrase)) {
                  phrase <- paste0(phrase, key_c, "=", category)
                } else {
                  phrase <- paste0(phrase, ", ", key_c, "=", category)
                }

                for (m in c:ncol(preds_modified)) {
                  if (m == ncol(preds_modified) & names(preds_modified)[m] %in% dict_key_c$match) {
                    next_index <- ncol(preds_modified) + 1
                  } else {
                    if (names(preds_modified)[m] %in% dict_key_c$match) {
                      next
                    } else {
                      next_index <- m
                      break
                    }
                  }
                }
              } else {
                if (is.null(phrase)) {
                  phrase <- paste0(phrase, match, "=", preds_modified[i, c])
                } else {
                  phrase <- paste0(phrase, ", ", match, "=", preds_modified[i, c])
                }
              }
            }
          }

          strata_preds[i] <- phrase
        }
      }

      preds$strata <- strata_preds
    }

    tibble_out$.pred[[k]] <- preds
  }

  return(tibble_out)
}

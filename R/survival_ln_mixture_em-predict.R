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
  
  tibble_out <- tibble::tibble(.pred = out)
  
  form <- model$blueprint$formula
  form <- stats::as.formula(paste0(paste(form)[[2]], " ~ ", paste(form)[[3]]))
  vars <- strsplit(paste(form)[[3]], " \\+ ")[[1]]
  
  if(!model$blueprint$intercept) {
    vars <- vars[vars != '0']
  }
  
  for(k in 1:length(tibble_out$.pred)) {
    preds <- tibble_out$.pred[[k]]
    
    if(all(vars != 'NULL')) {
      data <- model$data
      
      dict <- tibble::tibble(
        key = NA,
        match = NA,
        category = NA
      )
      
      predictors_data <- data |>
        dplyr::select(dplyr::all_of(vars))
      
      if(length(vars) == 1) {
        predictors_check_levels <- predictors_data
        names(predictors_check_levels) <- 'x'
        if(length(levels(droplevels(predictors_check_levels$x))) == 1) {        
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
        if(model$blueprint$intercept) {
          phrase <- NULL
          next_index <- 0
          for (c in 4:ncol(preds)) {
            if (c < next_index) {
              next
            }
            
            match <- names(preds)[c]
            key_c <- dict$key[dict$match == match]
            if (match %in% dict$match) {
              dict_key <- dict |>
                dplyr::filter(key == key_c)
              
              if (nrow(dict_key) > 2) { # categorical covariate with more than 2 levels
                category <- dict_key$category[!(dict_key$match %in% names(preds))]
                
                for (j in c:ncol(preds)) {
                  match_j <- names(preds)[j]
                  key_j <- dict$key[dict$match == match_j]
                  
                  if (key_j == key_c) {
                    if (preds[i, j] == 1) {
                      category <- dict_key$category[dict_key$match == match_j]
                    }
                  } else {
                    next_index <- j
                    break
                  }
                }
              } else {
                if (preds[i, c] == 1) {
                  category <- dict_key$category[dict_key$match == match]
                } else {
                  category <- dict_key$category[dict_key$match != match]
                }
              }
              
              phrase <- paste0(phrase, key_c, "=", category)
              
              if (c < ncol(preds)) {
                phrase <- paste0(phrase, ", ")
              }
            } else {
              phrase <- paste0(phrase, match, "=", preds[i, c])
              
              if (c < ncol(preds)) {
                phrase <- paste0(phrase, ", ")
              }
            }
          }
          strata_preds[i] <- phrase
        } else {
          phrase <- NULL
          next_index <- 0
          if(length(unique(dict$key)) == 1) {
            cols_remove <- dict |> 
              dplyr::group_by(key) |> 
              dplyr::slice(-n()) |> 
              dplyr::pull(match)
          } else {
            cols_remove <- NULL
            for(ky in unique(dict$key)) {
              if(sum(startsWith(names(preds), ky)) == 
                 nrow(dict |> dplyr::filter(key == ky))) {
                cols_remove <- c(cols_remove, dict$match[dict$key == ky][1])
              }
            }
          }
          
          preds_modified <- preds |>
            dplyr::select(-dplyr::all_of(cols_remove))
          
          for (c in 3:ncol(preds_modified)) {
            if (c < next_index) {
              next
            }
            
            match <- names(preds_modified)[c]
            key_c <- dict$key[dict$match == match]
            if (match %in% dict$match) {
              dict_key <- dict |>
                dplyr::filter(key == key_c)
              
              if (nrow(dict_key) > 2) { # categorical covariate with more than 2 levels
                category <- dict_key$category[!(dict_key$match %in% names(preds_modified))]
                
                for (j in c:ncol(preds_modified)) {
                  match_j <- names(preds_modified)[j]
                  key_j <- dict$key[dict$match == match_j]
                  
                  if (key_j == key_c) {
                    if (preds_modified[i, j] == 1) {
                      category <- dict_key$category[dict_key$match == match_j]
                    }
                  } else {
                    next_index <- j
                    break
                  }
                }
              } else {
                if (preds_modified[i, c] == 1) {
                  category <- dict_key$category[dict_key$match == match]
                } else {
                  category <- dict_key$category[dict_key$match != match]
                }
              }
              
              phrase <- paste0(phrase, key_c, "=", category)
              
              if (c < ncol(preds_modified)) {
                phrase <- paste0(phrase, ", ")
              }
            } else {
              phrase <- paste0(phrase, match, "=", preds_modified[i, c])
              
              if (c < ncol(preds_modified)) {
                phrase <- paste0(phrase, ", ")
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
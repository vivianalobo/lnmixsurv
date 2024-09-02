#' Function used to calculate the binary cross entropy (also known as log loss), given by
#' H_p(q) = -1/n * sum(S * log(\hat{S}) + (1 - S) * log(1 - \hat{S})), where S is the empirical survival and \hat{S} is the fitted survival.
#' 
#' `binary_cross_entropy()` is used to calculate the binary cross entropy of a predictions object, preferably the object `$preds` returned from `plot_fit_on_data()`.
#'
#' @param preds The `$preds` object from `plot_fit_on_data()` applied to the model. If not this one, should be a prediction `tibble` with the columns `time`, `strata` (if applicable), `estimate`, `.pred_survival`, `n.risk` (also `chain`, if necessary). It's important that the quantities `estimate` and `.pred_survival` are calculated for the same `time` and `strata`. It's highly recommended to simply use the object `$preds` returned from the function `plot_fit_on_data()`.
#'
#' 
#' @param nobs The number of observations used to fit the model. Can be ignored if `threshold` is set to 0. To easily calculate this value, use the function `nobs()` applied to the model object.
#'
#' @param threshold Numeric value between 0 and 1. Times with `n.risk` below threshold * nobs will be ignored. Default is 0.005 (0.5%).
#' 
#' @returns A `tibble` with the following columns:
#' - `strata`: The stratas used to fit the model.
#' - `mean_entropy`: Mean binary cross entropy across times.
#' - `chain`: The chain of the Bayesian model (only if necessary).
#'
#' @export
binary_cross_entropy <- function(preds, nobs = NULL, threshold = 0.005) {
  # Checks
  if (!is_tibble(preds)) {
    stop("The input 'preds' must be a tibble.")
  }
  
  if (!is.numeric(threshold) || threshold < 0 || threshold > 1) {
    stop("The input 'threshold' must be a numeric value between 0 and 1.")
  }
  
  if (!is.null(nobs) && (!is.numeric(nobs) || nobs <= 0)) {
    stop("The input 'nobs' must be a positive numeric value.")
  }
  
  if (is.null(nobs) && threshold != 0) {
    stop("The input 'nobs' must be provided if 'threshold' is different from 0.")
  }
  
  # Defining the numeric threshold (nt)
  if (threshold == 0) {
    nt <- 0
  } else {
    nt <- threshold * nobs
  }
  
  # Calculating the binary cross entropy
  out <- tibble::tibble(
    time = rep(NA, nrow(preds)),
    strata = rep(NA, nrow(preds)),
    entropy = rep(NA, nrow(preds)),
    chain = rep(NA, nrow(preds))
  )
  
  for(i in 1:nrow(preds)) {
    out$time[i] <- preds$time[i]
    
    if ('strata' %in% names(preds)) {
      out$strata[i] <- preds$strata[i]
    }
    
    if ('chain' %in% names(preds)) {
      out$chain[i] <- preds$chain[i]
    }
    
    if(preds$n.risk[i] > nt) {
      out$entropy[i] <- single_entropy(preds$estimate[i], preds$.pred_survival[i])
    }
  }
  
  if ('strata' %in% names(preds)) {
    if ('chain' %in% names(preds)) {
      out <- out |> 
        dplyr::group_by(strata, chain) |> 
        dplyr::summarise(mean_entropy = mean(entropy, na.rm = TRUE))
    } else {
      out <- out |> 
        dplyr::group_by(strata) |> 
        dplyr::summarise(mean_entropy = mean(entropy, na.rm = TRUE))
    }
  } else {
    if ('chain' %in% names(preds)) {
      out <- out |> 
        dplyr::group_by(chain) |> 
        dplyr::summarise(mean_entropy = mean(entropy, na.rm = TRUE))
    } else {
      out <- out |> 
        dplyr::summarise(mean_entropy = mean(entropy, na.rm = TRUE))
    }
  }
  
  return(out)
}

single_entropy <- function(x, y) {
  -(x * log(y) + (1 - x) * log(1 - y))
}

#' Function used to calculate some distance metrics between the predicted survival and the observed survival.
#' `fit_metrics()` is used to calculate distance metrics between empirical and fitted survival for a predictions object, preferably the object `$preds` returned from `plot_fit_on_data()`.
#'
#' @param preds The `$preds` object from `plot_fit_on_data()` applied to the model. If not this one, should be a prediction `tibble` with the columns `time`, `strata` (if applicable), `estimate`, `.pred_survival`, `n.risk` (also `chain`, if necessary). It's important that the quantities `estimate` and `.pred_survival` are calculated for the same `time` and `strata`. It's highly recommended to simply use the object `$preds` returned from the function `plot_fit_on_data()`.
#'
#' @param nobs The number of observations used to fit the model. Can be ignored if `threshold` is set to 0. To easily calculate this value, use the function `nobs()` applied to the model object.
#'
#' @param threshold Numeric value between 0 and 1. Times with `n.risk` below threshold * nobs will be ignored. Default is 0.005 (0.5%). Important because the distance metrics may be too big if calculated in intervals without sufficient observations to be estimated.
#' 
#' @returns A `tibble` with the following columns:
#' - `strata`: The stratas used to fit the model (if necessary).
#' - `n_strata`: The number of observations in the strata.
#' - `chain`: The chain of the Bayesian model (only if necessary).
#' - `metric`: Which metric is being calculated.
#' - `value`: Value for the metric.
#' 
#' For now, the following metrics are available and will be included:
#' - `MSE`: Mean Squared Error (the less the better).
#' - `MAE`: Mean Absolute Error (the less the better).
#' - `Hellinger Distance`: Hellinger distance, sometimes called Jeffreys distance (the less the better).
#' - `KS Distance`: Kolmogorov-Smirnov distance (the less the better).
#'
#' @export
fit_metrics <- function(preds, nobs = NULL, threshold = 0.005) {
  # Checks
  if (!tibble::is_tibble(preds)) {
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
  
  n.risk <- metric <- value <- nvalid <- value_mean <- value_max <- NULL
  
  # Defining the numeric threshold (nt)
  if (threshold == 0) {
    nt <- 0
  } else {
    nt <- threshold * nobs
  }
  
  epsilon <- 1e-10
  
  out <- list()
  
  for(i in 1:nrow(preds)) {
    time <- preds$time[i]
    
    if ('strata' %in% names(preds)) {
      strata <- preds$strata[i]
      n_strata <- preds |> 
        dplyr::filter(strata == preds$strata[i]) |> 
        dplyr::pull(n.risk) |>
        max()
    } else {
      strata <- NA
      n_strata <- preds |> 
        dplyr::pull(n.risk) |>
        max()
    }
    
    if ('chain' %in% names(preds)) {
      chain <- preds$chain[i]
    } else {
      chain <- NA
    }
    
    if(preds$n.risk[i] > nt) {
      out_tibble <- tibble::tibble(time = rep(time, 4),
                                   strata = rep(strata, 4),
                                   n_strata = rep(n_strata, 4),
                                   chain = rep(chain, 4),
                                   metric = c('MSE', 'MAE', 'Hellinger Distance', 'KS Distance'),
                                   value = c((preds$estimate[i] - preds$.pred_survival[i])^2,
                                             abs(preds$estimate[i] - preds$.pred_survival[i]),
                                             (sqrt(preds$estimate[i]) - sqrt(preds$.pred_survival[i]))^2,
                                             abs(preds$estimate[i] - preds$.pred_survival[i])))
      out[[i]] <- out_tibble
    }
  }
  
  out <- dplyr::bind_rows(out)
  
  if ('strata' %in% names(preds)) {
    if ('chain' %in% names(preds)) {
      out <- out |> 
        stats::na.omit() |> 
        dplyr::group_by(strata, chain, metric) |> 
        dplyr::summarise(value_mean = mean(value, na.rm = TRUE),
                         value_max = max(value, na.rm = TRUE),
                         nvalid = dplyr::n(),
                         n_strata = mean(n_strata)) |> 
        dplyr::mutate(value = ifelse(metric == 'Hellinger Distance',
                                     (1/sqrt(2)) * sqrt(nvalid * value_mean),
                                     ifelse(metric == 'KS Distance',
                                            value_max,
                                            value_mean))) |> 
        dplyr::select(-nvalid, -value_max, -value_mean)
    } else {
      out <- out |> 
        dplyr::select(-chain) |> 
        stats::na.omit() |> 
        dplyr::group_by(strata, metric) |> 
        dplyr::summarise(value_mean = mean(value, na.rm = TRUE),
                         value_max = max(value, na.rm = TRUE),
                         nvalid = dplyr::n(),
                         n_strata = mean(n_strata)) |> 
        dplyr::mutate(value = ifelse(metric == 'Hellinger Distance',
                                     (1/sqrt(2)) * sqrt(nvalid * value_mean),
                                     ifelse(metric == 'KS Distance',
                                            value_max,
                                            value_mean))) |> 
        dplyr::select(-nvalid, -value_max, -value_mean)
    }
  } else {
    if ('chain' %in% names(preds)) {
      out <- out |> 
        dplyr::select(-strata) |>
        stats::na.omit() |>
        dplyr::group_by(chain, metric) |> 
        dplyr::summarise(value_mean = mean(value, na.rm = TRUE),
                         value_max = max(value, na.rm = TRUE),
                         nvalid = dplyr::n(),
                         n_strata = mean(n_strata)) |> 
        dplyr::mutate(value = ifelse(metric == 'Hellinger Distance',
                                     (1/sqrt(2)) * sqrt(nvalid * value_mean),
                                     ifelse(metric == 'KS Distance',
                                            value_max,
                                            value_mean))) |> 
        dplyr::select(-nvalid, -value_max, -value_mean)
    } else {
      out <- out |> 
        dplyr::select(-strata, -chain) |> 
        stats::na.omit() |>
        dplyr::group_by(metric) |> 
        dplyr::summarise(value_mean = mean(value, na.rm = TRUE),
                         value_max = max(value, na.rm = TRUE),
                         nvalid = dplyr::n(),
                         n_strata = mean(n_strata)) |> 
        dplyr::mutate(value = ifelse(metric == 'Hellinger Distance',
                                     (1/sqrt(2)) * sqrt(nvalid * value_mean),
                                     ifelse(metric == 'KS Distance',
                                            value_max,
                                            value_mean))) |> 
        dplyr::select(-nvalid, -value_max, -value_mean)
    }
  }
  
  return(out)
}

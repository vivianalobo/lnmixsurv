#' Function used to quick visualize the fitted values (survival estimate) on the data used to fit the model.
#'
#' `plot_fit_on_data_em()` estimates survival/hazard for the data the model was fitted on and plots the results.
#'
#' @param model A `survival_ln_mixture_em` object.
#'
#' @param type A character string specifying the type of plot. The default is "survival", but can be "hazard".
#'
#' @export
plot_fit_on_data_em <- function(model, type = "survival") {
  # Checks
  if (!inherits(model, "survival_ln_mixture_em") &
      !inherits(model, "survival_ln_mixture")) {
    stop("The object should be a survival_ln_mixture_em or a survival_ln_mixture object.")
  }
  
  if(inherits(model, "survival_ln_mixture_em") & 
     inherits(model, "survival_ln_mixture")) {
    stop("The object should just one of two: survival_ln_mixture_em or survival_ln_mixture. Can't be both at the sime time.")
  }
  
  if (!type %in% c("survival", "hazard")) {
    stop("The type should be either 'survival' or 'hazard'.")
  }
  
  # -----
  # Predictions
  
  # which is the type of model fitted
  if(inherits(model, "survival_ln_mixture_em")) {
    data <- model$data
    form <- model$blueprint$formula
    form <- stats::as.formula(paste0(paste(form)[[2]], " ~ ", paste(form)[[3]]))
    vars <- strsplit(paste(form)[[3]], " \\+ ")[[1]]
    if(!model$blueprint$intercept) {
      vars <- vars[vars != '0']
    }
    
    time <- NA
    estimate <- NA
    .eval_time <- NA
    .pred_survival <- NA
    .pred_hazard <- NA
    
    km <- survival::survfit(form, data) |>
      broom::tidy()
    
    eval_time <- seq(0, max(km$time))
    
    if(all(vars != 'NULL')) {
      preds <- stats::predict(model,
                              dplyr::distinct(data |>
                                                dplyr::select(dplyr::all_of(vars))),
                              type = type,
                              eval_time = eval_time
      ) |>
        tidyr::unnest(.pred)
    } else {
      preds <- stats::predict(model,
                              data.frame(val = NA),
                              type = type,
                              eval_time = eval_time
      ) |>
        tidyr::unnest(.pred)
    }
    
    
    if(type == 'survival') {
      if(all(vars != 'NULL')) {
        gg <- ggplot() +
          geom_step(aes(x = time, y = estimate, color = strata),
                    data = km, alpha = 0.5) +
          geom_line(aes(x = .eval_time, y = .pred_survival, color = strata),
                    data = preds) +
          theme_bw() +
          labs(x = 'Time', y = 'Survival')
      } else {
        gg <- ggplot() +
          geom_step(aes(x = time, y = estimate),
                    data = km, alpha = 0.5) +
          geom_line(aes(x = .eval_time, y = .pred_survival),
                    data = preds) +
          theme_bw() +
          labs(x = 'Time', y = 'Survival')
      }
    } else {
      # yet to be implemented hazard
    }
  } else {
    # yet to be implemented bayesian
  }
  
  return(list(preds = preds,
              ggplot = gg))
}

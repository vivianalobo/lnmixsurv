#' Function used to quick visualize the fitted values (survival estimate) on the data used to fit the model.
#'
#' `plot_fit_on_data_em()` estimates survival/hazard for the data the model was fitted on and plots the results.
#'
#' @param model A `survival_ln_mixture_em` object.
#'
#' @param type A character string specifying the type of plot. The default is "survival", but can be "hazard".
#'
#' @export
plot_fit_on_data_em <- function(model, type = 'survival') {
  # check the object class
  if (!inherits(model, "survival_ln_mixture_em")) {
    stop("The object should be a survival_ln_mixture_em object.")
  }
  
  # check type
  if (!type %in% c('survival', 'hazard')) {
    stop("The type should be either 'survival' or 'hazard'.")
  }
  
  form <- model$blueprint$formula
  form <- stats::as.formula(paste0(paste(form)[[2]], ' ~ ', paste(form)[[3]]))
  vars <- strsplit(paste(form)[[3]], ' \\+ ')[[1]]
  
  data <- model$data
  
  dict <- tibble::tibble(key = NA,
                         match = NA,
                         category = NA)
  
  predictors <- data |> 
    dplyr::select(dplyr::all_of(vars))
  
  time <- NA
  estimate <- NA
  .eval_time <- NA
  .pred_survival <- NA
  .pred_hazard <- NA
  
  for(c in 1:ncol(predictors)) {
    if(is.factor(predictors |> 
                 dplyr::select(dplyr::all_of(c)) |> 
                 dplyr::pull())) {
      
      key <- names(predictors)[c]
      var <- predictors |> 
        dplyr::select(dplyr::all_of(c)) |> 
        dplyr::pull()
      
      for(i in levels(var)) {
        match <- paste0(key, i)
        category <- i
        dict <- dplyr::bind_rows(dict, tibble::tibble(key = key,
                                                      match = match,
                                                      category = category))
      } 
    }
  }
  
  dict <- dict |> 
    dplyr::slice(-1)
  
  km <- survival::survfit(form, data) |> 
    broom::tidy()
  
  eval_time <- seq(0, max(km$time))
  
  preds <- stats::predict(model,
                          dplyr::distinct(data |> 
                                            dplyr::select(dplyr::all_of(vars))),
                          type = type,
                          eval_time = eval_time) |> 
    tidyr::unnest(.pred)
  
  strata_preds <- rep(NA, nrow(preds))
  
  for(i in 1:nrow(preds)) {
    phrase <- NULL
    next_index <- 0
    for(c in 4:ncol(preds)) {
      if(c < next_index) {
        next
      }
      
      match <- names(preds)[c]
      key_c <- dict$key[dict$match == match]
      if(match %in% dict$match) {
        dict_key <- dict |> 
          dplyr::filter(key == key_c)
        
        if(nrow(dict_key) > 2) { # categorical covariate with more than 2 levels
          category <- dict_key$category[!(dict_key$match %in% names(preds))]
          
          for(j in c:ncol(preds)) {
            match_j <- names(preds)[j]
            key_j <- dict$key[dict$match == match_j]
            
            if(key_j == key_c) {
              if(preds[i, j] == 1) {
                category <- dict_key$category[dict_key$match == match_j]
              }
            } else {
              next_index <- j
              break
            }
          }
        } else {
          if(preds[i, c] == 1) {
            category <- dict_key$category[dict_key$match == match]
          } else {
            category <- dict_key$category[dict_key$match != match]
          }
        }
        
        phrase <- paste0(phrase, key_c, '=', category)
        
        if(c < ncol(preds)) {
          phrase <- paste0(phrase, ', ')
        }
      } else {
        phrase <- paste0(phrase, match, '=', preds[i, c])
        
        if(c < ncol(preds)) {
          phrase <- paste0(phrase, ', ')
        }
      }
    }
    strata_preds[i] <- phrase
  }
  
  preds$strata <- strata_preds
  
  if(type == 'survival') {
    gg <- ggplot() +
      geom_step(aes(x = time, y = estimate, color = strata), data = km,
                alpha = 0.5) +
      geom_line(aes(x = .eval_time, y = .pred_survival, color = strata), data = preds) +
      theme_bw() +
      labs(x = 'Time',
           y = 'Survival')
  } else {
    gg <- ggplot() +
      geom_step(aes(x = time, y = estimate, color = strata), data = km,
                alpha = 0.5) +
      geom_line(aes(x = .eval_time, y = .pred_hazard, color = strata), data = preds) +
      theme_bw() +
      labs(x = 'Time',
           y = 'Survival')
  }
  
  return(gg)
}

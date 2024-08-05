#' Function used to quick visualize the fitted values (survival estimate) on the data used to fit the model (via EM algorithm or Gibbs).
#'
#' `plot_fit_on_data()` estimates survival/hazard for the data the model was fitted on and plots the results.
#'
#' @param model A `survival_ln_mixture` or `survival_ln_mixture_em` object.
#'
#' @param data A `data.frame()` or `tibble()` containing the data used to fit the model.
#'
#' @param type A character string specifying the type of plot. The default is "survival", but can be "hazard".
#'
#' @param interval A character string specifying the type of interval to be plotted. The default is "none", but can be "credible". The EM algorithm does not provide confidence intervals and this parameter is only support for the Bayesian version (`survival_ln_mixture` object).
#'
#' @param level A numeric value between 0 and 1 specifying the level of the confidence interval. The default is 0.95.
#'
#' @export
plot_fit_on_data <- function(model, data, type = "survival", interval = "none",
                             level = 0.95) {
  # Checks
  if (!inherits(model, "survival_ln_mixture_em") &
    !inherits(model, "survival_ln_mixture")) {
    stop("The object should be a survival_ln_mixture_em or a survival_ln_mixture object.")
  }

  if (inherits(model, "survival_ln_mixture_em") &
    inherits(model, "survival_ln_mixture")) {
    stop("The object should just one of two: survival_ln_mixture_em or survival_ln_mixture. Can't be both at the sime time.")
  }

  if (!type %in% c("survival", "hazard")) {
    stop("The type should be either 'survival' or 'hazard'.")
  }

  # Check level numeric
  if (!is.numeric(level)) {
    stop("The level should be a numeric value.")
  }

  # Check level
  if (level < 0 | level > 1) {
    stop("The level should be between 0 and 1.")
  }

  # -----
  # Predictions
  # -----

  form <- model$blueprint$formula
  form <- stats::as.formula(paste0(paste(form)[[2]], " ~ ", paste(form)[[3]]))
  vars <- strsplit(paste(form)[[3]], " \\+ ")[[1]]
  if (!model$blueprint$intercept) {
    vars <- vars[vars != "0"]
  }

  time <- estimate <- .eval_time <- .pred_survival <- .pred_hazard <- .pred_upper <- .pred_lower <- credible_ribbon <- facet_chain <- line_layer <- step_layer <- guides_gg <- NA

  km <- survival::survfit(form, data) |>
    broom::tidy()

  eval_time <- seq(0, max(km$time))

  # which is the type of model fitted
  if (inherits(model, "survival_ln_mixture_em")) {
    # warn that the interval and level
    if (interval != "none") {
      warning("The EM algorithm does not provide confidence intervals. It will be ignored.")
    }

    if (all(vars != "NULL")) {
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

    if (type == "survival") {
      if (all(vars != "NULL")) {
        gg <- ggplot() +
          geom_step(aes(x = time, y = estimate, color = strata),
            data = km, alpha = 0.5
          ) +
          geom_line(aes(x = .eval_time, y = .pred_survival, color = strata),
            data = preds
          ) +
          theme_bw() +
          labs(x = "Time", y = "Survival")
      } else {
        gg <- ggplot() +
          geom_step(aes(x = time, y = estimate),
            data = km, alpha = 0.5
          ) +
          geom_line(aes(x = .eval_time, y = .pred_survival),
            data = preds
          ) +
          theme_bw() +
          labs(x = "Time", y = "Survival")
      }
    } else {
      km <- join_empirical_hazard(km)

      if (all(vars != "NULL")) {
        gg <- ggplot() +
          geom_line(aes(x = time, y = hazard_estimate, color = strata),
            data = km, alpha = 0.5
          ) +
          geom_line(aes(x = .eval_time, y = .pred_hazard, color = strata),
            data = preds
          ) +
          theme_bw() +
          labs(x = "Time", y = "Hazard")
      } else {
        gg <- ggplot() +
          geom_line(aes(x = time, y = hazard_estimate),
            data = km, alpha = 0.5
          ) +
          geom_line(aes(x = .eval_time, y = .pred_hazard),
            data = preds
          ) +
          theme_bw() +
          labs(x = "Time", y = "Hazard")
      }
    }
  } else {
    if (all(vars != "NULL")) {
      preds <- stats::predict(model,
        dplyr::distinct(data |>
          dplyr::select(dplyr::all_of(vars))),
        type = type,
        eval_time = eval_time,
        interval = interval,
        level = level
      ) |>
        tidyr::unnest(.pred)
    } else {
      preds <- stats::predict(model,
        data.frame(val = NA),
        type = type,
        eval_time = eval_time,
        interval = interval,
        level = level
      ) |>
        tidyr::unnest(.pred)
    }

    if (type == "survival") {
      labs_gg <- labs(x = "Time", y = "Survival")
      if (posterior::nchains(model$posterior) > 1) {
        facet_chain <- facet_wrap(~chain)
        guides_gg <- guides(linetype = "none")
        if (all(vars != "NULL")) {
          step_layer <- geom_step(aes(x = time, y = estimate, color = strata), data = km, alpha = 0.5)
          line_layer <- geom_line(aes(x = .eval_time, y = .pred_survival, color = strata, linetype = chain), data = preds)
          if (interval == "credible") {
            credible_ribbon <- geom_ribbon(aes(x = .eval_time, ymin = .pred_lower, ymax = .pred_upper, fill = strata, linetype = chain), data = preds, alpha = 0.3)
          } else {
            credible_ribbon <- NULL
          }
        } else {
          step_layer <- geom_step(aes(x = time, y = estimate), data = km, alpha = 0.5)
          line_layer <- geom_line(aes(x = .eval_time, y = .pred_survival, linetype = chain), data = preds)

          if (interval == "credible") {
            credible_ribbon <- geom_ribbon(aes(x = .eval_time, ymin = .pred_lower, ymax = .pred_upper, linetype = chain), data = preds, alpha = 0.3)
          } else {
            credible_ribbon <- NULL
          }
        }
      } else {
        guides_gg <- NULL
        facet_chain <- NULL
        if (all(vars != "NULL")) {
          line_layer <- geom_line(aes(x = .eval_time, y = .pred_survival, color = strata), data = preds)
          step_layer <- geom_step(aes(x = time, y = estimate, color = strata), data = km, alpha = 0.5)

          if (interval == "credible") {
            credible_ribbon <- geom_ribbon(aes(x = .eval_time, ymin = .pred_lower, ymax = .pred_upper, fill = strata), data = preds, alpha = 0.3)
          } else {
            credible_ribbon <- NULL
          }
        } else {
          line_layer <- geom_line(aes(x = .eval_time, y = .pred_survival), data = preds)
          step_layer <- geom_step(aes(x = time, y = estimate), data = km, alpha = 0.5)

          if (interval == "credible") {
            credible_ribbon <- geom_ribbon(aes(x = .eval_time, ymin = .pred_lower, ymax = .pred_upper), data = preds, alpha = 0.3)
          } else {
            credible_ribbon <- NULL
          }
        }
      }

      gg <- ggplot() +
        step_layer +
        line_layer +
        credible_ribbon +
        labs_gg +
        theme_bw() +
        guides_gg +
        facet_chain
    } else {
      km <- join_empirical_hazard(km)
      labs_gg <- labs(x = "Time", y = "Hazard")

      if (posterior::nchains(model$posterior) > 1) {
        facet_chain <- facet_wrap(~chain)
        guides_gg <- guides(linetype = "none")
        if (all(vars != "NULL")) {
          step_layer <- geom_line(aes(x = time, y = hazard_estimate, color = strata), data = km, alpha = 0.5)
          line_layer <- geom_line(aes(x = .eval_time, y = .pred_hazard, color = strata, linetype = chain), data = preds)

          if (interval == "credible") {
            credible_ribbon <- geom_ribbon(aes(x = .eval_time, ymin = .pred_lower, ymax = .pred_upper, fill = strata, linetype = chain), data = preds, alpha = 0.3)
          } else {
            credible_ribbon <- NULL
          }
        } else {
          step_layer <- geom_line(aes(x = time, y = hazard_estimate), data = km, alpha = 0.5)
          line_layer <- geom_line(aes(x = .eval_time, y = .pred_hazard, linetype = chain), data = preds)

          if (interval == "credible") {
            credible_ribbon <- geom_ribbon(aes(x = .eval_time, ymin = .pred_lower, ymax = .pred_upper, linetype = chain), data = preds, alpha = 0.3)
          } else {
            credible_ribbon <- NULL
          }
        }
      } else {
        guides_gg <- NULL
        facet_chain <- NULL
        if (all(vars != "NULL")) {
          line_layer <- geom_line(aes(x = .eval_time, y = .pred_hazard, color = strata), data = preds)
          step_layer <- geom_line(aes(x = time, y = hazard_estimate, color = strata), data = km, alpha = 0.5)

          if (interval == "credible") {
            credible_ribbon <- geom_ribbon(aes(x = .eval_time, ymin = .pred_lower, ymax = .pred_upper, fill = strata), data = preds, alpha = 0.3)
          } else {
            credible_ribbon <- NULL
          }
        } else {
          line_layer <- geom_line(aes(x = .eval_time, y = .pred_hazard), data = preds)
          step_layer <- geom_line(aes(x = time, y = hazard_estimate), data = km, alpha = 0.5)

          if (interval == "credible") {
            credible_ribbon <- geom_ribbon(aes(x = .eval_time, ymin = .pred_lower, ymax = .pred_upper), data = preds, alpha = 0.3)
          } else {
            credible_ribbon <- NULL
          }
        }
      }

      gg <- ggplot() +
        step_layer +
        line_layer +
        credible_ribbon +
        labs_gg +
        theme_bw() +
        guides_gg +
        facet_chain
    }
  }

  return(list(
    preds = preds,
    ggplot = gg
  ))
}

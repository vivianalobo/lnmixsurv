make_survival_reg_survival_ln_mixture <- function() {
  parsnip::set_model_engine(
    "survival_reg",
    mode = "censored regression",
    eng = "survival_ln_mixture"
  )
  parsnip::set_dependency(
      "survival_reg",
      eng = "survival_ln_mixture", pkg = "persistencia"
  )

  parsnip::set_fit(
    model = "survival_reg",
    eng = "survival_ln_mixture",
    mode = "censored regression",
    value = list(
      interface = "formula",
      protect = c("formula", "data"),
      func = c(pkg = "persistencia", fun = "survival_ln_mixture"),
      defaults = list()
    )
  )

  parsnip::set_encoding(
    model = "survival_reg",
    eng = "survival_ln_mixture",
    mode = "censored regression",
    options = list(
      predictor_indicators = "traditional",
      compute_intercept = FALSE,
      remove_intercept = FALSE,
      allow_sparse_x = FALSE
    )
  )

  parsnip::set_pred(
    model = "survival_reg",
    eng = "survival_ln_mixture",
    mode = "censored regression",
    type = "survival",
    value = list(
      pre = NULL,
      post = NULL,
      func = c(fun = "predict"),
      args =
        list(
          object = quote(object$fit),
          new_data = quote(new_data),
          type = "survival",
          time = quote(time)
        )
    )
  )

  parsnip::set_pred(
    model = "survival_reg",
    eng = "survival_ln_mixture",
    mode = "censored regression",
    type = "hazard",
    value = list(
      pre = NULL,
      post = NULL,
      func = c(fun = "predict"),
      args =
        list(
          object = quote(object$fit),
          new_data = quote(new_data),
          type = "hazard",
          time = quote(time)
        )
    )
  )
}
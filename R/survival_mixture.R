survival_mixture <- function(mode = "censored regression") {
  # Check for correct mode
  if (mode != "censored regression") {
    rlang::abort("`mode` should be 'censored regression'")
  }

  # Save some empty slots for future parts of the specification
  parsnip::new_model_spec(
    "survival_mixture",
    args = NULL,
    eng_args = NULL,
    mode = mode,
    method = NULL,
    engine = NULL
  )
}

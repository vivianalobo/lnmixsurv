#' @export
print.survival_ln_mixture <- function(model, digits = NULL, ...) {
  if (is.null(digits)) digits <- getOption("digits")
  fixed <- tidy.survival_ln_mixture(model, effects = "fixed")
  auxiliary <- tidy.survival_ln_mixture(model, effects = "auxiliary")
  cat("survival_ln_mixture")
  cat("\n formula:     ", extract_formula(model))
  cat("\n observations:", nobs(model))
  cat("\n predictors:  ", npredictors(model))
  cat("\n------\n")
  print(data.frame(fixed[, -1], row.names = fixed$term), digits = digits)
  cat("\nAuxiliary parameter(s):\n")
  print(data.frame(auxiliary[, -1], row.names = auxiliary$term), digits = digits)
}

#' Funcao auxiliar para calcular intervalo de credibilidade usando quantis.
#' @noRd
interval <- function(x, conf.level) {
  ret <- unname(quantile(x, probs = c(1 - conf.level, conf.level)))
  return(c("conf.low" = ret[1], "conf.high" = ret[2]))
}

#' @export
tidy.survival_ln_mixture <- function(x, # nolint: object_name_linter.
                                     effects = "fixed",
                                     conf.int = FALSE, # nolint: object_name_linter.
                                     conf.level = 0.9, # nolint: object_name_linter.
                                     digits = NULL,
                                     ...) {
  rlang::arg_match(effects, c("fixed", "auxiliary"))

  vars <- c()
  if ("fixed" %in% effects) {
    vars <- c(vars, "beta_a", "beta_b")
  }
  if ("auxiliary" %in% effects) {
    vars <- c(vars, "phi_a", "phi_b", "theta_a")
  }
  measures <- c("estimate" = median, "std.error" = mad)
  if (conf.int) {
    measures <- c(measures, c("interval" = function(x) interval(x, conf.level)))
  }
  post <- posterior::subset_draws(x$posterior, variable = vars)
  post <- posterior::merge_chains(post)
  ret <- posterior::summarise_draws(post, measures)
  names(ret)[1] <- "term"
  ret
}

#' Tidying method for a Lognormal Mixture model.
#'
#' These method tidy the estimates from `survival_ln_mixture` fits into a summary.
#'
#' @param x Fitted model object.
#' @param effects A character vector including one or more of `"fixed"` and `"auxiliary`.
#' @param conf.int If `TRUE` columns for lower (`conf.low`) and upper (`conf.high`) bounds
#' of the posterior uncertainty intervals are included.
#' @param conf.level A number between 0 and 1 indicating the desired probability mass to include in the
#' intervals. Only used if `conf.int = TRUE`.
#' @param digits How many significant digits should be displayed?
#' @param ... Not used.
#'
#' @return
#'
#' A `data.frame` without rownames. When `effects="fixed"` (the default), tidy.survival_ln_mixutre
#' returns one row for each coefficient for each component of the mixture with three columns:
#' \item{term}{The name of the corresponding term in the model.}
#' \item{estimate}{A point estimate of the coefficient (posterior median).}
#' \item{std.error}{A standard error for the point estimate based on
#' \code{\link[stats]{mad}}. See the \emph{Uncertainty estimates} section in
#' \code{\link[rstanarm]{print.stanreg}} for more details.}
#'
#' Setting \code{effects="auxiliary"} will select the remaining parameters:
#'
#' \item{phi_a}{Dispersion parameter for the first componente of the mixture.}
#' \item{phi_b}{Dispersion parameter for the second componente of the mixture.}
#' \item{theta_a}{The weigth of the first component of the mixture.}
#'
#' @examples
#'
#' require(survival)
#' lung$sex <- factor(lung$sex)
#' set.seed(1)
#' mod2 <- survival_ln_mixture(Surv(time, status == 2) ~ sex, lung, intercept = TRUE)
#' tidy(mod2)
#' tidy(mod2, conf.int = TRUE)
#' tidy(mod2, effects = c("fixed", "auxiliary"), conf.int = TRUE)
#'
#' @export
tidy.survival_ln_mixture <- function(x, # nolint: object_name_linter.
                                     effects = "fixed",
                                     conf.int = FALSE, # nolint: object_name_linter.
                                     conf.level = 0.9, # nolint: object_name_linter.
                                     digits = NULL,
                                     ...) {
  rlang::arg_match(effects, c("fixed", "auxiliary"))
  rlang::check_dots_empty(...)
  vars <- c()
  if ("fixed" %in% effects) {
    vars <- c(vars, "beta_a", "beta_b")
  }
  if ("auxiliary" %in% effects) {
    vars <- c(vars, "phi_a", "phi_b", "theta_a")
  }
  measures <- c("estimate" = stats::median, "std.error" = stats::mad)
  if (conf.int) {
    measures <- c(measures, c("interval" = function(x) interval(x, conf.level)))
  }
  post <- posterior::subset_draws(x$posterior, variable = vars)
  post <- posterior::merge_chains(post)
  ret <- posterior::summarise_draws(post, measures)
  names(ret)[1] <- "term"
  ret
}

#' Funcao auxiliar para calcular intervalo de credibilidade usando quantis.
#' @noRd
interval <- function(x, conf.level) { # nolint: object_name_linter.
  ret <- unname(stats::quantile(x, probs = c(1 - conf.level, conf.level)))
  return(c("conf.low" = ret[1], "conf.high" = ret[2]))
}

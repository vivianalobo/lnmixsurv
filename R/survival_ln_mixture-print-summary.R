# print.survival_ln_mixture <- function() {

# }

#' @export
summary.survival_ln_mixture <- function(object) {
    posterior::summarise_draws(object$posterior)
    invisible(object)
}

#' @importFrom generics tidy
#' @export
generics::tidy
#' @importFrom stats nobs
#' @export
stats::nobs

# nocov start

.onLoad <- function(libname, pkgname) {
    make_survival_reg_survival_ln_mixture()
}

# nocov end

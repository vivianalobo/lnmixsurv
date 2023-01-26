#' @importFrom generics tidy
#' @export
generics::tidy
#' @importFrom stats nobs
#' @export
stats::nobs

.onLoad <- function(libname, pkgname) {
  make_survival_reg_survival_ln_mixture()
}

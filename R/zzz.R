#' @importFrom generics glance tidy
#' @export
generics::glance
#' @export
generics::tidy


.onLoad <- function(libname, pkgname) {
  make_survival_reg_survival_ln_mixture()
}

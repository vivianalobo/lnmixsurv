#' @importFrom generics tidy
#' @export
generics::tidy

.onLoad <- function(libname, pkgname) {
  make_survival_reg_survival_ln_mixture()
}

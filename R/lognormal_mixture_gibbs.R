#' LOGNORMAL CENSORED GIBBS SCHEME
#'
#' Gibbs estimation of a lognormal mixture model.
#' (include here model equation)
#'
#' @param x design matrix
#' @param y responde vector in the original scale. all must be positive, since
#' the log transfomation
#' will be taken.
#' @param M number of iterations
#' @param cens binary vector of censored observations
#' (0 = censored ; 1 = observed)
#'
#' @return list with the posterior distribution of each parameter.
#'
#' @export
#'
lognormal_mixture_gibbs <- function(x, y, m, delta, intercepto = F, valor_inicial_beta = 0) {
    covs <- factor(x)
    cov_names <- levels(covs)
    if (intercepto) {
        x <- stats::model.matrix(~covs)
    } else {
        x <- stats::model.matrix(~ covs - 1)
    }

    t1 <- Sys.time()
    res <- lognormal_mixture_gibbs_cpp(y, x, delta, m, valor_inicial_beta)
    t2 <- Sys.time()
    tempo <- t2 - t1
    message(glue::glue("Rodado em {round(tempo,2)} {units(tempo)}."))
    res$betaA <- t(res$beta[, 1, ])
    res$betaB <- t(res$beta[, 2, ])
    res$beta <- NULL
    colnames(res$betaA) <- cov_names
    colnames(res$betaB) <- cov_names


    return(lapply(res, posterior::as_draws_matrix))
}

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
lognormal_mixture_gibbs <- function(x, y, delta, intercepto = FALSE, numero_iteracoes = 1000, numero_cadeias = 1, valor_inicial_beta = 0) { # nolint: line_length_linter.
    ncores <- parallel::detectCores(logical = FALSE)
    if (numero_cadeias > ncores) {
        # TODO: incluir acentos em ascii
        msg <- paste0(
            "O numero de cadeias ", numero_cadeias,
            " deve ser menor ou igual ao o numero de cores disponiveis",
            ncores, "."
        )
        stop(msg)
    }
    covs <- factor(x)
    cov_names <- levels(covs)
    if (intercepto) {
        x <- stats::model.matrix(~covs)
    } else {
        x <- stats::model.matrix(~ covs - 1)
    }

    t1 <- Sys.time()
    if (numero_cadeias == 1) {
        res <- lognormal_mixture_gibbs_cpp(y, x, delta, numero_iteracoes, valor_inicial_beta)
    } else {
        res <- parallel_lognormal_mixture_gibbs_cpp(y, x, delta, numero_iteracoes, numero_cadeias, valor_inicial_beta)
    }

    t2 <- Sys.time()
    tempo <- t2 - t1
    message(glue::glue("Rodado em {round(tempo,2)} {units(tempo)}."))
    names(res) <- c("betaA", "betaB", "phiA", "phiB", "theta")
    dimnames(posteriori$betaA) <- list(NULL, NULL, cov_names)
    dimnames(posteriori$betaB) <- list(NULL, NULL, cov_names)


    return(lapply(res, posterior::as_draws_array))
}

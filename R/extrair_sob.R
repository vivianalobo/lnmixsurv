sob_lognormal <- function(t, m, sigma) {
  s <- stats::pnorm((-log(t) + m) / sigma)
  return(s)
}

sob_lognormal_mix <- function(t, ma, mb, sigmaa, sigmab, theta) {
  s <- sob_lognormal(t, ma, sigmaa) * theta + sob_lognormal(t, mb, sigmab) * (1 - theta)
  return(s)
}


falha_lognormal_mix <- function(t, ma, mb, sigmaa, sigmab, theta, sob_mix) {
  sob_mix <- sob_lognormal_mix(t, ma, mb, sigmaa, sigmab, theta)
  dlnorm_mix <- theta * stats::dlnorm(t, ma, sigmaa) + (1 - theta) * stats::dlnorm(t, mb, sigmab)
  return(dlnorm_mix / sob_mix)
}


falha_lognormal_mix_eficiente <- function(t, ma, mb, sigmaa, sigmab, theta, sob_mix) {
  dlnorm_mix <- theta * stats::dlnorm(t, ma, sigmaa) + (1 - theta) * stats::dlnorm(t, mb, sigmab)
  return(dlnorm_mix / sob_mix)
}

funcao_sobrevivencia_risco_cenario <- function(post, t, cenario) {
  qntd_iteracoes <- nrow(post[[1]])
  cenarios <- colnames(post[[1]])

  if (!cenario %in% cenarios) stop("Cenario invalido.")

  sob <- matrix(NA, nrow = qntd_iteracoes, ncol = length(t))
  risco <- matrix(NA, nrow = qntd_iteracoes, ncol = length(t))
  betaA <- posterior::subset_draws(post$betaA, variable = cenario)
  betaB <- posterior::subset_draws(post$betaB, variable = cenario)
  for (iteration in seq(qntd_iteracoes)) {
    ma <- betaA[iteration]
    mb <- betaB[iteration]
    sigmaa <- sqrt(1 / post$phiA[iteration])
    sigmab <- sqrt(1 / post$phiB[iteration])
    theta <- post$theta[iteration]
    sob[iteration, ] <- sob_lognormal_mix(t, ma, mb, sigmaa, sigmab, theta)
    risco[iteration, ] <- falha_lognormal_mix_eficiente(
      t, ma, mb, sigmaa, sigmab, theta, sob[iteration, ]
    )
  }
  sob <- apply(sob, 2, stats::median)
  risco <- apply(risco, 2, stats::median)

  return(data.table::data.table(time = t, surv = sob, wx = risco))
}

#' @export
funcao_sobrevivencia_risco <- function(post, t) {
  cenarios <- colnames(post[[1]])
  all_sob <- list()

  for (cenario in cenarios) {
    all_sob[[cenario]] <- funcao_sobrevivencia_risco_cenario(post, t, cenario)
  }
  return(data.table::rbindlist(all_sob, idcol = "cenarios"))
}

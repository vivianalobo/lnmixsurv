sob_lognormal <- function(t, m, sigma) {
    s <- stats::pnorm((-log(t) + m) / sigma)
    return(s)
}

sob_lognormal_mix <- function(t, ma, mb, sigmaa, sigmab, theta) {
    s <- sob_lognormal(t, ma, sigmaa) * theta + sob_lognormal(t, mb, sigmab) * (1 - theta)
    return(s)
}


falha_lognormal_mix <- function(t, ma, mb, sigmaa, sigmab, theta) {
    sob_mix <- sob_lognormal_mix(t, ma, mb, sigmaa, sigmab, theta)
    dlnorm_mix <- theta * stats::dlnorm(t, ma, sigmaa) + (1 - theta) * stats::dlnorm(t, mb, sigmab)
    return(dlnorm_mix / sob_mix)
}


falha_lognormal_mix_eficiente <- function(t, ma, mb, sigmaa, sigmab, theta, sob_mix) {
    dlnorm_mix <- theta * stats::dlnorm(t, ma, sigmaa) + (1 - theta) * stats::dlnorm(t, mb, sigmab)
    return(dlnorm_mix / sob_mix)
}

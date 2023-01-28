require(truncnorm)
gerar_dados <- function(n = 1000, percentual_censura = 0.4, theta = 0.6) {
    beta <- rbind(c(3.3, 0.5), c(4, .8)) # coeficientes de regressão para as componentes
    V <- c(0.3, 0.0392207) # sigma^2
    sigma <- sqrt(V) # sigma

    x <- rbinom(n, 1, 0.4) # covariavel (dicotomica)
    I <- rbinom(n, 1, theta) # indicadora de primeira componente de mistura
    delta <- rbinom(n, 1, 1 - percentual_censura) # status (censura ou falha)

    mu <- cbind(1, x) %*% t(beta) # parametro mu por componente

    y1 <- rnorm(n, mu[, 1], sigma[1])
    y2 <- rnorm(n, mu[, 2], sigma[2])
    y <- I * y1 + (1 - I) * y2

    for (i in which(delta == 0)) { # truncando e censurado
        y[i] <- I[i] * rtruncnorm(1, a = -Inf, b = y[i], mean = mu[i, 1], sd = sigma[1]) +
            (1 - I[i]) * rtruncnorm(1, a = -Inf, b = y[i], mean = mu[i, 2], sd = sigma[2])
    }
    true_vals <- c(
        "beta_a[1]" = beta[1, 1],
        "beta_a[2]" = beta[1, 2],
        "beta_b[1]" = beta[2, 1],
        "beta_b[2]" = beta[2, 2],
        "phi_a" = 1 / V[1],
        "phi_b" = 1 / V[2],
        "theta_a" = theta
    )
    return(list(data = data.frame(y = exp(y), x = factor(x), delta = delta), true_vals = true_vals))
}

set.seed(1)
sim_data <- gerar_dados(10000)

usethis::use_data(sim_data, overwrite = TRUE)

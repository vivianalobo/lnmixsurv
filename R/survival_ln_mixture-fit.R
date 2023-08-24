#' Lognormal mixture model
#'
#' `survival_ln_mixture()` fits a lognormal mixture model, as described in (referencia artigo viviana).
#' Colocar mais detalhes sobre o modelo (como a equacao) em details.
#'
#' @param formula A formula specifying the outcome terms on the left-hand side,
#' and the predictor terms on the right-hand side. The outcome must be a [survival::Surv]
#' object.
#'
#' @param data A __data frame__ containing both the predictors and the outcome.
#'
#' @param intercept A logical. Should an intercept be included in the processed data?
#'
#' @param iter A positive integer specifying the number of iterations for each chain (including warmup).
#' 
#' @param em_iter A positive integer specifying the number of iterations for the EM algorithm. The EM algorithm is performed before the Gibbs sampler to find better initial values for the chains. On simulations, values lower than 200 seems to work nice.
#'
#' @param warmup A positive integer specifying the number of warmup (aka burnin) iterations per chain.
#' The number of warmup iterations should be smaller than iter.
#'
#' @param thin A positive integer specifying the period for saving samples.
#'
#' @param chains A positive integer specifying the number of Markov chains.
#'
#' @param cores Ignored. Parallel runs are disabled.
#' 
#' @param numero_componentes number of mixture componentes >= 2.
#'
#' @param proposal_variance The value used at the distribution for e0, hyperparameter of the Dirichlet prior, has the form of Gamma(proposal_variance, proposal_variance*G). It affects how distant the proposal values will be from the actual value. Large values of the proposal_variance may be problematic, since the hyperparameter e0 is sampled using a Metropolis-Hasting algorithm and may take long to converge. The code is implemented so the initial value of proposal_variance does not affect the convergence too much, since it's changed through the iterations to sintonize the variance, ensuring an acceptance ratio of proposal values between 17% and 25%, which seems to be optimal on our tests.
#' 
#' @param show_progress Indicates if the code shows the progress of the EM algorithm and the Gibbs Sampler.
#' 
#' @param ... Not currently used, but required for extensibility.
#'
#' @note Categorical predictors must be converted to factors before the fit,
#' otherwise the predictions will fail.
#'
#' @return
#'
#' A `survival_ln_mixture` object, which is a list with the following componentes:
#'
#' \item{posterior}{A [posterior::draws_matrix] with the posterior of the parameters of the model.}
#' \item{nobs}{A integer holding the number of observations used to generate the fit.}
#' \item{blueprint}{The blueprint component of the output of [hardhat::mold]}
#'
#'
#' @examples
#'
#' # Formula interface
#' library(survival)
#' set.seed(1)
#' mod <- survival_ln_mixture(Surv(time, status == 2) ~ NULL, lung, intercept = TRUE)
#'
#' @export
survival_ln_mixture <- function(formula, data, intercept = TRUE, iter = 1000, warmup = floor(iter / 10), thin = 1, chains = 1, cores = 1, numero_componentes = 2, proposal_variance = 2, show_progress = FALSE, em_iter = 150, ...) {
  rlang::check_dots_empty(...)
  UseMethod("survival_ln_mixture")
}

#' @export
#' @rdname survival_ln_mixture
survival_ln_mixture.default <- function(formula, ...) {
  stop("`survival_ln_mixture()` is not defined for a '", class(formula)[1], "'.", call. = FALSE)
}

# Formula method

#' @export
#' @rdname survival_ln_mixture
survival_ln_mixture.formula <- function(formula, data, intercept = TRUE, ...) {
  blueprint <- hardhat::default_formula_blueprint(intercept = intercept)
  processed <- hardhat::mold(formula, data, blueprint = blueprint)
  survival_ln_mixture_bridge(processed, ...)
}

# ------------------------------------------------------------------------------
# Bridge

survival_ln_mixture_bridge <- function(processed, ...) {
  predictors <- as.matrix(processed$predictors)
  outcome <- processed$outcome[[1]]

  if (!survival::is.Surv(outcome)) {
    rlang::abort("Response must be a survival object (created with survival::Surv)")
  }
  if (attr(outcome, "type") != "right") rlang::abort("Only right-censored data allowed")

  outcome_times <- outcome[, 1]
  outcome_status <- outcome[, 2]

  fit <- survival_ln_mixture_impl(predictors, outcome_times, outcome_status, ...)

  new_survival_ln_mixture(
    posterior = fit$posterior,
    nobs = fit$nobs,
    predictors_name = fit$predictors_name,
    mixture_groups = fit$mixture_groups,
    blueprint = processed$blueprint
  )
}


# ------------------------------------------------------------------------------
# Implementation

survival_ln_mixture_impl <- function(predictors, outcome_times, 
                                     outcome_status, iter = 1000,
                                     warmup = floor(iter / 10), thin = 1,
                                     chains = 1, cores = 1, 
                                     numero_componentes = 2,
                                     proposal_variance = 1,
                                     show_progress = FALSE,
                                     em_iter = 150) {
  
  number_of_predictors <- ncol(predictors)
  
  if (number_of_predictors < 1) {
    rlang::abort(
      c(
        "The model must contain at least one predictor.",
        i = "When using outcome ~ NULL, intercept must be explicitly set to TRUE."
      )
    )
  }
  
  if (any(outcome_times == 0)) {
    rlang::abort("One or more events happened at time zero.")
  }
  
  if (cores != 1) warning("Argumento cores ignorado, rodando cadeias sequencialmente.")
  
  posterior_dist <- sequential_lognormal_mixture_gibbs(
    iter, em_iter, numero_componentes, chains, outcome_times, outcome_status,
    predictors, proposal_variance, show_progress)
  
  posterior_dist <- give_colnames(posterior_dist, colnames(predictors),
                                  numero_componentes)
  
  posterior_dist <- label_switch_one_chain(posterior_dist)
  
  posterior_dist <- permute_columns(posterior_dist)
  remover_menor_theta <- -which(
    colnames(posterior_dist) == colnames(
      posterior_dist |> 
        dplyr::select(dplyr::starts_with('eta_')))[numero_componentes])
  posterior_dist <- posterior_dist[, remover_menor_theta]
  
  posterior_dist <- posterior::as_draws_matrix(posterior_dist)
  posterior_dist <- posterior::subset_draws(posterior_dist, iteration = seq(from = warmup + 1, to = iter))
  posterior_dist <- posterior::thin_draws(posterior_dist, thin = thin)

  list(posterior = posterior_dist, 
       nobs = length(outcome_times),
       predictors_name = colnames(predictors), 
       mixture_groups = seq_len(numero_componentes))
}

#' Converte um field of cubes (lista de arrays) para uma array representando a posteriori de uma cadeia.
#' 
#' @param field field of cubes, exportado do rcpp. No R, é uma lista de arrays. A primeira array da lista possui
#' dimensao numero_covariaveis x numero_componetes x numero_iteracoes, e é referente ao parametro beta.
#' A segunda e terceira array da lista possuem dimensao 1 x numero_componetnes x numero_iteracoes e são referentes
#' aos parametros phi e theta respectivamente.
#' @param pred_names nome das variáveis preditoras. Deve ser um vetor com tamanho numero_covariaveis.
#' @param grupos nome dado as componetes (usualmente a, b, c)
#' 
#' @return Uma array de dimensão numero_iteracoes x 1 x numero_componetes * (2 + numero_covariaveis). Representa a 
#' distribuicao a posteriori de uma cadeia, num formato aceito pelo posterior::as_draws_matrix.
#' 
#' @noRd
field_of_cube2array <- function(field, pred_names, grupos){
  posterior_dist = abind::abind(field, along = 1)
  
  names_beta <- glue::glue_data(expand.grid(pred_names, grupos), "{Var1}_{Var2}")
  names_phi <- glue::glue("phi_{grupos}")
  names_theta <- glue::glue("theta_{grupos}")
  
  dimnames(posterior_dist)[[1]] <- c(pred_names, "phi", "theta")
  dimnames(posterior_dist)[[2]] <- grupos
  
  colnames = glue::glue_data(expand.grid(dimnames(posterior_dist)[1:2]), "{Var1}_{Var2}")
  
  posterior_dist = t(apply(posterior_dist, 3, c))
  posterior_dist = array(
    posterior_dist, 
    dim = c(nrow(posterior_dist), ncol(posterior_dist), 1), 
    dimnames = list(NULL, colnames, NULL)
  )
  posterior_dist = aperm(posterior_dist, c(1, 3, 2))[,,c(names_beta, names_phi, names_theta), drop = FALSE]
  
  return(posterior_dist)
}

#' Aplica fild_of_cubes2array em todas as cadeias do field_of_fields exportado pelo Rcpp.
#' 
#' @param field_of_field liste de lista de arrays exportada pelo Rcpp. Cada elemento da lista é uma lista contendo 
#' 3 arrays de dimensoes numero_covariaveis x numero_componentes x numero_iteracoes, 1 x numero_componetes x numero_iteracoes
#' e 1 x numero_componentes x numero_iteracoes representando os parametros beta, phi e theta respectivamente.
#' @param pred_names nome das variáveis preditoras. Deve ser um vetor com tamanho numero_covariaveis.
#' @param grupos nome dado as componetes (usualmente a, b, c)
#' 
#' @return objecto da classe draws_matrix.
#' 
#' @noRd
field_of_fied2array <- function(field_of_field, pred_names, grupos){
  list_of_arrays = lapply(field_of_field, function(x) field_of_cube2array(x, pred_names, grupos))
  arrays_joined = abind::abind(list_of_arrays, along = 2) 
  return(arrays_joined)
}

#' corrige o problema do label switch para uma cadeia da posteriori, ordenando grupos por proporções de mistura
#' 
#' @param posterior_dist uma matriz de dimensao numero_iteracoes x numero_componentes de amostras a posteriori

#' @return matriz de dimensão numero_iteracoes x numero_componetes mas com os labels
#' reorganizados de forma que os etas das componetes são ordenados de forma decrescente.
#' 
#' @noRd
label_switch_one_chain <- function(posterior_dist){
  obj <- dplyr::as_tibble(posterior_dist)
  
  obj_etas <- obj |>
    dplyr::select(dplyr::starts_with('eta_'))
  
  etas_median <- as.numeric(apply(obj_etas, 2, stats::median))
  
  etas_order <- order(etas_median, decreasing = T)
  
  new_obj <- NULL
  
  for(j in 1:length(etas_order)) {
    sub_obj <- obj |>
      dplyr::select(dplyr::ends_with(as.character(etas_order[j])))
    
    for(c in 1:ncol(sub_obj)) {
      names(sub_obj)[c] <- paste0(
        substr(names(sub_obj)[c],
               1, nchar(names(sub_obj)[c]) - 1), j)
    }
    
    new_obj <- dplyr::bind_cols(new_obj, sub_obj)
  }
  
  return(new_obj)
}

#' corrige o problema do label switch para todas as cadeias da posteriori.
#' 
#' @param posterior_dist uma array de dimensao numero_iteracoes x numero_cadeias x numero_componetes * (2 + numero_covariaveis)
#' @param pred_names nome das variaveis preditoras
#' @param grupos nome dado as componetes (usualmente a, b, c)
#' 
#' @return uma array de dimensao numero_iteracoes x numero_cadeias x numero_componetes * (2 + numero_covariaveis) mas com os labels
#' reorganizados de forma que os thetas das componetes são ordenados de forma decrescente. 
#' 
#' @noRd
label_switch_all <- function(posterior_dist, pred_names, grupos){
  for(i in seq_len(dim(posterior_dist)[2])){
    posterior_dist[,i,] = label_switch_one_chain(posterior_dist[,i,], pred_names, grupos)
  }
  return(posterior_dist)
}

#' Nomeia as colunas da amostra a posteriori
#' 
#' @param posterior distribuição a posteriori amostrada
#' 
#' @param predictors_names nome das variáveis preditoras. Deve ser um vetor com tamanho numero_covariaveis.
#' @param numero_componentes número de componentes envolvidos no ajuste
#' 
#' @return matriz
#' 
#' @noRd
give_colnames <- function(posterior, predictors_names, numero_componentes) {
  number_params <- length(predictors_names)
  new_names <- NULL
  
  for(i in 1:numero_componentes) {
    for(j in 1:3) { # de 1 a 3 porque 3 grupos de parâmetros são ajustados no modelo: uma proporção de mistura, uma precisão e um grupo de efeitos das covariáveis
      if(j == 1) { # primeiro grupo, eta: efeitos das covariáveis
        for(c in 1:number_params) { 
          new_names <- c(new_names,
                         paste0(predictors_names[c], '_', i))
        }
      } else if (j == 2) { # segundo grupo, phi: precisão
        new_names <- c(new_names,
                       paste0('phi_', i))
      } else { # terceiro grupo, j = 3, proporções de mistura
        new_names <- c(new_names,
                       paste0('eta_', i))
      }
    }
  }
  
  posterior_dist <- posterior
  colnames(posterior_dist) <- new_names
  
  return(posterior_dist)
}

#' Permuta as colunas para ficar de acordo: primeiro efeitos das covariáveis dos grupos, depois precisões e, por fim, etas
#' 
#' @param posterior distribuição a posteriori amostrada
#' 
#' @return matriz
#' 
#' @noRd
permute_columns <- function(posterior) {
  posterior_dist <- dplyr::as_tibble(posterior)
  posterior_dist <- dplyr::bind_cols(posterior |> 
                                       dplyr::select(
                                         -tidyselect::starts_with('eta'),
                                         -tidyselect::starts_with('phi')),
                                     posterior |> 
                                       dplyr::select(
                                         tidyselect::starts_with('phi'),
                                         tidyselect::starts_with('eta')))
  
  return(posterior_dist)
}

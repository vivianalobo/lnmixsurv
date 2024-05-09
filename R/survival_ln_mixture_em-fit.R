#' Lognormal mixture model - Expectation-Maximization Algorithm
#'
#' `survival_ln_mixture_em()` fits an EM algorithm, as described in LOBO, Viviana GR; FONSECA, Thaís CO; ALVES, Mariane B. Lapse risk modeling in insurance: a Bayesian mixture approach. Annals of Actuarial Science, v. 18, n. 1, p. 126-151, 2024, for modelling mixtures of lognormal distributions applied to survival data.
#'
#' @param formula A formula specifying the outcome terms on the left-hand side,
#' and the predictor terms on the right-hand side. The outcome must be a [survival::Surv]
#' object.
#'
#' @param data A __data frame__ containing both the predictors and the outcome.
#'
#' @param iter A positive integer specifying the number of iterations for the EM algorithm.
#'
#' @param mixture_components number of mixture componentes >= 2.
#'
#' @param starting_seed Starting seed for the algorithm. If not specified by the user, uses a random integer between 1 and 2^28 This way we ensure, when the user sets a seed in R, that this is passed into the C++ code.
#'
#' @param intercept A logical. Should an intercept be included in the processed data?
#'
#' @param sparse Useful if the design matrix is sparse (most cases with categorical only regressors). Can save a lot of memory, allowing for huge data to be fitted.
#'
#' @param ... Not currently used, but required for extensibility.
#'
#' @export
survival_ln_mixture_em <- function(formula, data, intercept = TRUE, iter = 50, mixture_components = 2, starting_seed = sample(1:2^28, 1), sparse = FALSE, ...) {
  rlang::check_dots_empty(...)
  UseMethod("survival_ln_mixture_em")
}

#' @export
#' @rdname survival_ln_mixture_em
survival_ln_mixture_em.default <- function(formula, ...) {
  stop("`survival_ln_mixture_em()` is not defined for a '",
    class(formula)[1], "'.",
    call. = FALSE
  )
}

# Formula method
#' @export
#' @rdname survival_ln_mixture_em
survival_ln_mixture_em.formula <- function(formula, data, intercept = TRUE, ...) {
  blueprint <- hardhat::default_formula_blueprint(intercept = intercept)
  processed <- hardhat::mold(formula, data, blueprint = blueprint)
  survival_ln_mixture_em_bridge(processed, ...)
}

# ---------------- BRIDGE ----------------
survival_ln_mixture_em_bridge <- function(processed, ...) {
  # Verifications
  outcome <- processed$outcome[[1]]
  if (!survival::is.Surv(outcome)) {
    rlang::abort("Response must be a survival object (created with survival::Surv)")
  }

  if (attr(outcome, "type") != "right") {
    rlang::abort("Only right-censored data allowed")
  }

  predictors <- as.matrix(processed$predictors)
  outcome_times <- outcome[, 1]
  outcome_status <- outcome[, 2]

  fit <- survival_ln_mixture_em_impl(
    outcome_times, outcome_status,
    predictors, ...
  )

  new_survival_ln_mixture_em(
    em_iterations = fit$em_iterations,
    nobs = fit$nobs,
    predictors_name = fit$predictors_name,
    mixture_groups = fit$mixture_groups,
    blueprint = processed$blueprint
  )
}

# ---------------- IMPLEMENTATION ----------------
survival_ln_mixture_em_impl <- function(outcome_times, outcome_status,
                                        predictors, iter = 50,
                                        mixture_components = 2,
                                        starting_seed = sample(1:2^28, 1),
                                        sparse = FALSE) {
  # Verifications
  if (any(is.na(predictors))) {
    "There is one or more NA values in the predictors variable."
  }

  number_predictors <- ncol(predictors)

  if (number_predictors < 1) {
    rlang::abort(
      c("The model must contain at least one predictor.",
        i = "When using outcome ~ NULL, intercept must be explicitly set to TRUE."
      )
    )
  }

  if (any(outcome_times == 0)) {
    rlang::abort("One or more events happened at time zero.")
  }

  if (any(is.na(outcome_times))) {
    rlang::abort("There is one or more NA values at event times.")
  }

  if (any(is.na(outcome_status))) {
    rlang::abort("There is one or more NA values at the status")
  }

  if (iter <= 0 | (iter %% 1) != 0) {
    rlang::abort("The parameter iter should be a positive integer.")
  }

  if (starting_seed < 1 | starting_seed > 2^28 |
    (starting_seed %% 1) != 0) {
    rlang::abort("The starting seed should be a natural number between 1 and 2^28")
  }

  if (mixture_components <= 0 | (mixture_components %% 1) != 0) {
    rlang::abort("The parameter mixture_components should be a positive integer.")
  }

  # These next two lines seems to be unecessary but they are essencial to ensure
  # the reproducibility of the EM iterations on the Gibbs sampler. For an user,
  # interested only in using the EM, this is irrelevant.
  set.seed(starting_seed)

  seed <- sample(1:2^28, 1)

  matrix_em_iter <- lognormal_mixture_em_implementation(
    iter, mixture_components, outcome_times,
    outcome_status, predictors, seed, sparse
  )

  new_names <- NULL

  for (g in 1:mixture_components) {
    for (j in 1:3) {
      if (j == 1) {
        new_names <- c(
          new_names,
          paste0("eta_", g)
        )
      } else if (j == 2) {
        for (k in 0:(number_predictors - 1)) {
          new_names <- c(
            new_names,
            paste0("beta", k, "_", g)
          )
        }
      } else {
        new_names <- c(
          new_names,
          paste0("phi_", g)
        )
      }
    }
  }

  colnames(matrix_em_iter) <- new_names
  matrix_em_iter <- dplyr::bind_cols(
    matrix_em_iter,
    tibble::tibble(iter = 1:iter)
  )

  list(
    em_iterations = matrix_em_iter,
    number_iterations = iter,
    nobs = length(outcome_times),
    mixture_groups = mixture_components,
    predictors_name = colnames(predictors)
  )
}

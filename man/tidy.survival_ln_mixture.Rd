% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidiers.R
\name{tidy.survival_ln_mixture}
\alias{tidy.survival_ln_mixture}
\title{Tidying method for a Lognormal Mixture model.}
\usage{
\method{tidy}{survival_ln_mixture}(
  x,
  effects = "fixed",
  conf.int = FALSE,
  conf.level = 0.9,
  digits = NULL,
  ...
)
}
\arguments{
\item{x}{Fitted model object.}

\item{effects}{A character vector including one or more of \code{"fixed"} and \verb{"auxiliary}.}

\item{conf.int}{If \code{TRUE} columns for lower (\code{conf.low}) and upper (\code{conf.high}) bounds
of the posterior uncertainty intervals are included.}

\item{conf.level}{A number between 0 and 1 indicating the desired probability mass to include in the
intervals. Only used if \code{conf.int = TRUE}.}

\item{digits}{How many significant digits should be displayed?}

\item{...}{Not used.}
}
\value{
A \code{data.frame} without rownames. When \code{effects="fixed"} (the default), tidy.survival_ln_mixutre
returns one row for each coefficient for each component of the mixture with three columns:
\item{term}{The name of the corresponding term in the model.}
\item{estimate}{A point estimate of the coefficient (posterior median).}
\item{std.error}{A standard error for the point estimate based on
\code{\link[stats]{mad}}. See the \emph{Uncertainty estimates} section in
\code{\link[rstanarm]{print.stanreg}} for more details.}

Setting \code{effects="auxiliary"} will select the remaining parameters:

\item{phi_a}{Dispersion parameter for the first componente of the mixture.}
\item{phi_b}{Dispersion parameter for the second componente of the mixture.}
\item{theta_a}{The weigth of the first component of the mixture.}
}
\description{
These method tidy the estimates from \code{survival_ln_mixture} fits into a summary.
}
\examples{

require(survival)
lung$sex <- factor(lung$sex)
set.seed(1)
mod2 <- survival_ln_mixture(Surv(time, status == 2) ~ sex, lung, intercept = TRUE)
tidy(mod2)
tidy(mod2, conf.int = TRUE)
tidy(mod2, effects = c("fixed", "auxiliary"), conf.int = TRUE)

}

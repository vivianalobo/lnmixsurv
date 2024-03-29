% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{sim_data}
\alias{sim_data}
\title{Simulated lognormal mixture data.}
\format{
\subsection{\code{sim_data}}{

A list with two componentes:
\itemize{
\item $data: A data frame with 10,000 rows and 3 columns:
\describe{
\item{y}{observed survival time}
\item{delta}{event indicator. 1 == event, 0 == censored.}
\item{x}{binary covariate}
}
\item $true_vals: A named vector with the true values used to generate the data.
}
}
}
\usage{
sim_data
}
\description{
A simulated dataset with 10000 observations from a lognormal mixutre
model with 2 componentes.
}
\keyword{datasets}

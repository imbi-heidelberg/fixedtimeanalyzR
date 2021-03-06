% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exp_surv-data.R
\docType{data}
\name{exp_surv}
\alias{exp_surv}
\title{Exponential survival data frame}
\format{
A dataframe with variables \code{time}, \code{group} and \code{status}.
}
\source{
Simulation by Lukas D. Sauer.
}
\usage{
data(exp_surv)
}
\description{
Simulated survival data with exponential event times and exponential
 censoring times.
}
\examples{
library(ggplot2)
library(ggfortify)
library(survival)

data(exp_surv)
exp_fit <- survfit(Surv(time, status) ~ group, data=exp_surv)
autoplot(exp_fit)
}
\keyword{datasets}

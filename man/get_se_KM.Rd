% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_utils.R
\name{get_se_KM}
\alias{get_se_KM}
\title{Standard error of the Kaplan-Meier survival}
\usage{
get_se_KM(data = NULL, se_KM = NULL, t = NULL)
}
\arguments{
\item{data}{either null object or a \code{\link[survival]{survfit}} object or a data frame that has
the variables \code{time} specifying time to event or time to censoring,
\code{status} specifying the censoring status and \code{group} specifying the
 different groups to be compared.}

\item{se_KM}{either a null object or a numeric vector containing standard
errors for one or two groups}

\item{t}{the point in time at which the survival rate will be calculated}
}
\value{
A numeric vector containing the survival rates for the different
groups at time t.
}
\description{
This function returns the standard error of the Kaplan-Meier estimate of the
survival rate at a fixed point in time for two groups. There are three
possibilities for data input, (1) data with two groups and a
null variable se_KM, (2) data with one group and a numeric vector
se_KM of length 1, (3) a null variable data and a numeric vector se_KM of
length 2. In case (1), it is a wrapper function of the standard error from the \code{\link{survival}}
package. In case (2), it also uses this function, and then concatenates it
with the numeric vector se_KM. In case (3), it simply returns the numeric
vector se_KM.
}
\examples{
data(exp_surv)
get_se_KM(exp_surv, t=1)

library(dplyr)
exp_surv \%>\% filter(group == 1) -> exp_surv_grp1
get_se_KM(exp_surv_grp1, se_KM=0.11, t=1)

get_se_KM(se_KM=c(0.11,0.12))

}

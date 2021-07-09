#' Exponential survival data.
#'
#' Simulated survival data with exponential event times and exponential
#'  censoring times.
#'
#' @docType data
#'
#' @usage data(exp_surv)
#'
#' @format A dataframe with variables \code{time}, \code{group} and \code{status}.
#'
#' @keywords datasets
#'
#'
#' @source Simulation by Lukas D. Sauer.
#'
#' @examples
#' library(ggplot2)
#' library(survival)
#'
#' data(exp_surv)
#' exp_fit <- survfit(Surv(time, status) ~ group, data=exp_surv)
#' autoplot(exp_fit) + stat_function(fun=function(x){1-pexp(x, log(1/p1))}
#'
#' \donttest{iplotCurves(phe, times)}
"exp_surv"

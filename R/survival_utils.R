# ------------------------------------------------------------------------------
# UTILITY FUNCTIONS FOR SURVIVAL ANALYSIS

#' Kaplan-Meier estimate of the survival rate.
#'
#' This function returns the Kaplan-Meier estimate of the survival rate at a
#' fixed point in time. In fact, it is a wrapper function of the survival
#' estimate from the
#' \code{\link[https://www.rdocumentation.org/packages/survival/versions/2.11-4/topics/survfit]{survival}}
#'  package.
#'
#' @param data a dataset that has the variables \code{time} specifying time to
#' event or time to censoring, \code{status} specifying the censoring status and
#' \code{group} specifying the different groups to be compared.
#' @param t the point in time at which the survival rate will be calculated
#'
#' @return A numeric vector containing the survival rates for the different
#' groups at time t.
#'
#' @examples
#' data(exp_surv)
#' get_surv_KM(exp_surv, 1)
#'
#' @export
get_surv_KM <- function(data, t){
  # TODO: Ist das überhaupt in Ordnung, dass ich hier extend = TRUE gewählt habe?
  surv_KM = summary(survfit(formula = Surv(time, status) ~ group, data=sample), times=t, extend = TRUE)$surv
  return(surv_KM)
}

#-------------------------------------------------------------------------------
# Return the Kaplan-Meier estimate of the standard error of the
# survival rate at a fixed time.
#
# sample -- a dataset with columns status, time and treatment
# t -- a number specifying the time at which survival rate shall be calculated
get_se_KM <- function(sample, t){
  se_KM = summary(survfit(formula = Surv(time, status) ~ treatment, data=sample), times=t, extend = TRUE)$std.err
  return(se_KM)
}

get_sigma_KM <- function(surv_KM, se_KM){
  return(se_KM/surv_KM)
}

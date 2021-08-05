# ------------------------------------------------------------------------------
# UTILITY FinUNCTIONS FOR SURVIVAL ANALYSIS

#' Kaplan-Meier estimate of the survival rate.
#'
#' This function returns the Kaplan-Meier estimate of the survival rate at a
#' fixed point in time. In fact, it is a wrapper function of the survival
#' estimate from the
#' \code{\link[https://www.rdocumentation.org/packages/survival/versions/2.11-4/topics/survfit]{survival}}
#'  package.
#'
#' @param data either a survfit object or a data frame that has the variables
#' \code{time} specifying time to event or time to censoring, \code{status}
#' specifying the censoring status and \code{group} specifying the different
#' groups to be compared.
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
  if(!is(data, "survfit")){
    data_survfit <- get_survfit(data)
    return(summary(data_survfit, times=t, extend = TRUE)$surv)
  }
  return(summary(data, times=t, extend = TRUE)$surv)
}

#' Standard error of the Kaplan-Meier survival rate.
#'
#' This function returns the standard error of the Kaplan-Meier estimate of the
#' survival rate at a fixed point in time. In fact, it is a wrapper function of
#' the survival estimate from the
#' \code{\link[https://www.rdocumentation.org/packages/survival/versions/2.11-4/topics/survfit]{survival}}
#'  package.
#'
#' @param data either a survfit object or a data frame that has the variables
#' \code{time} specifying time to event or time to censoring, \code{status}
#' specifying the censoring status and \code{group} specifying the different
#' groups to be compared.
#' @param t the point in time at which the standard error will be calculated
#'
#' @return A numeric vector containing the standard errors for the different
#' groups at time t.
#'
#' @examples
#' data(exp_surv)
#' get_se_KM(exp_surv, 1)
#'
#' @export
get_se_KM <- function(data, t){
  if(!is(data, "survfit")){
    data_survfit <- get_survfit(data)
    return(summary(data_survfit, times=t, extend = TRUE)$std.err)
  }
  return(summary(data, times=t, extend = TRUE)$std.err)
}


#' Square root of Greenwood's formula.
#'
#' This function returns the square root of the quotient of the Kaplan-Meier
#' estimate of the survival rate by its standard error. This quotient can be
#' calculated by a formula known as Greenwood's formula. This helper function is
#' simply introduced for better readability when comparing the implemented
#' functions to the ones from Klein et al.
#'
#' @param data either a survfit object or a data frame that has the variables
#' \code{time} specifying time to event or time to censoring, \code{status}
#' specifying the censoring status and \code{group} specifying the different
#' groups to be compared.
#' @param t the point in time at which Greenwood's formula will be calculated.
#'
#' @return A numeric vector containing the square root of Greenwood's formula for
#'  the different groups at time t.
#'
#' @examples
#' data(exp_surv)
#' get_sigma_KM(exp_surv, 1)
#'
#' @export
get_sigma_KM <- function(data, t){
  se_KM = get_se_KM(data, t)
  surv_KM = get_surv_KM(data, t)
  return(se_KM/surv_KM)
}

#' Function for generating survfit object.
#'
#' This function returns a survfit object that has been generated from a
#' survival data set.
#'
#' @param data a data frame that has the variables \code{time} specifying time to
#' event or time to censoring, \code{status} specifying the censoring status and
#' \code{group} specifying the different groups to be compared.
#'
#' @return A survfit object.
#'
#' @examples
#' data(exp_surv)
#' get_survfit(exp_surv)
#'
#' @export
get_survfit <- function(data){
  if(!("time" %in% names(data) && "status" %in% names(data)
       && "group" %in% names(data))) {
    stop("Data set must contain the variables time, status and group.")
  }
  # check for numeric time variable.
  if (!is.numeric(data$time)) {
    stop("Variable time must be a numeric vector.")
  }
  return(survfit(formula = Surv(time, status) ~ group, data=data))
}

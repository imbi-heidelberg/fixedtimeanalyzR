# ------------------------------------------------------------------------------
# UTILITY FUNCTIONS FOR SURVIVAL ANALYSIS

#' Kaplan-Meier estimate of the survival rate
#'
#' This function returns the Kaplan-Meier estimate of the survival rate at a
#' fixed point in time. In fact, it is a wrapper function of the survival
#' estimate from the \code{\link{survival}} package.
#'
#' @param data either a \code{\link[survival]{survfit}} object or a data frame that has the variables
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
#' get_surv_KM_from_data(exp_surv, 1)
#'
#' @export
get_surv_KM_from_data <- function(data, t){
  # TODO: Ist das überhaupt in Ordnung, dass ich hier extend = TRUE gewählt habe?
  if(!methods::is(data, "survfit")){
    data_survfit <- get_survfit(data)
    return(summary(data_survfit, times=t, extend = TRUE)$surv)
  }
  return(summary(data, times=t, extend = TRUE)$surv)
}

#' Standard error of the Kaplan-Meier survival rate
#'
#' This function returns the standard error of the Kaplan-Meier estimate of the
#' survival rate at a fixed point in time. In fact, it is a wrapper function of
#' the survival estimate from the \code{\link{survival}} package.
#'
#' @param data either a \code{\link[survival]{survfit}} object or a data frame that has the variables
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
#' get_se_KM_from_data(exp_surv, 1)
#'
#' @export
get_se_KM_from_data <- function(data, t){
  if(!methods::is(data, "survfit")){
    data_survfit <- get_survfit(data)
    return(summary(data_survfit, times=t, extend = TRUE)$std.err)
  }
  return(summary(data, times=t, extend = TRUE)$std.err)
}

#' Kaplan-Meier estimate of the survival rate
#'
#' This function returns the Kaplan-Meier estimate of the survival rate at a
#' fixed point in time for two groups. There are three possibilities for data
#' input, (1) data with two groups and a null variables surv_KM,
#' (2) data with one group and a numeric vector surv_KM of length 1,
#' (3) a null variable data and a numeric vector surv_KM of length 2. In case
#'  (1), it is a wrapper function of the survival estimate from the \code{\link{survival}}
#' package. In case (2), it also uses this function, and then concatenates it
#' with the numeric vector surv_KM. In case (3), it simply returns the numeric
#' vector surv_KM
#'
#' @param data either a null object or a \code{\link[survival]{survfit}} object or a data frame that has
#' the variables \code{time} specifying time to event or time to censoring,
#' \code{status} specifying the censoring status and \code{group} specifying the
#'  different groups to be compared.
#' @param surv_KM either a null object or a numeric vector containing survival
#' rates for one or two groups
#' @param t the point in time at which the survival rate will be calculated
#'
#' @return A numeric vector containing the survival rates for the different
#' groups at time t.
#'
#' @examples
#' data(exp_surv)
#' get_surv_KM(exp_surv, t=1)
#'
#' library(dplyr)
#' exp_surv %>% filter(group == 1) -> exp_surv_grp1
#' get_surv_KM(exp_surv_grp1, surv_KM=0.5, t=1)
#'
#' get_surv_KM(surv_KM=c(0.5,0.6))
#'
#' @export
get_surv_KM <- function(data = NULL, surv_KM = NULL, t = NULL){
  n_groups = 0
  if(methods::is(data, "data.frame")){
    n_groups = length(levels(as.factor(data$group)))
  }
  else if(methods::is(data, "survfit")){
    n_groups = length(data$n)
  }
  if(!is.null(data) & is.null(surv_KM) & n_groups == 2){
    return(get_surv_KM_from_data(data, t))
  }
  else if(length(surv_KM) == 1 & (n_groups == 1 | n_groups == 0) & is.numeric(surv_KM)){
    if(n_groups == 0){warning("No group variable was detected in the data frame. Assuming that there is only one group.")}
    return(c(get_surv_KM_from_data(data,t), surv_KM))
  } else if(length(surv_KM) == 2 & is.null(data) & is.numeric(surv_KM)){
    if(!is.null(t)){
      warning("If you want to compare two numeric values of survival
            the time variable t is not used.")
    }
    return(surv_KM)
  } else{
    stop("Invalid input. There are three possibilities for data input,
        (1) a data frame or a survfit object 'data' with two groups and null variables 'surv_KM' and
        'se_KM', (2) a data frame or a survfit object 'data' with one group and numeric vectors 'surv_KM' and
        'se_KM' of length 1, (3) a null variable 'data' and numeric vectors 'surv_KM'
        and 'se_KM' of length 2.")
  }
}
#' Standard error of the Kaplan-Meier survival
#'
#' This function returns the standard error of the Kaplan-Meier estimate of the
#' survival rate at a fixed point in time for two groups. There are three
#' possibilities for data input, (1) data with two groups and a
#' null variable se_KM, (2) data with one group and a numeric vector
#' se_KM of length 1, (3) a null variable data and a numeric vector se_KM of
#' length 2. In case (1), it is a wrapper function of the standard error from the \code{\link{survival}}
#' package. In case (2), it also uses this function, and then concatenates it
#' with the numeric vector se_KM. In case (3), it simply returns the numeric
#' vector se_KM.
#'
#' @param data either null object or a \code{\link[survival]{survfit}} object or a data frame that has
#' the variables \code{time} specifying time to event or time to censoring,
#' \code{status} specifying the censoring status and \code{group} specifying the
#'  different groups to be compared.
#' @param se_KM either a null object or a numeric vector containing standard
#' errors for one or two groups
#' @param t the point in time at which the survival rate will be calculated
#'
#' @return A numeric vector containing the survival rates for the different
#' groups at time t.
#'
#' @examples
#' data(exp_surv)
#' get_se_KM(exp_surv, t=1)
#'
#' library(dplyr)
#' exp_surv %>% filter(group == 1) -> exp_surv_grp1
#' get_se_KM(exp_surv_grp1, se_KM=0.11, t=1)
#'
#' get_se_KM(se_KM=c(0.11,0.12))
#'
#' @export
get_se_KM <- function(data = NULL, se_KM = NULL, t = NULL){
  n_groups = 0
  if(methods::is(data, "data.frame")){
    n_groups = length(levels(as.factor(data$group)))
  }
  else if(methods::is(data, "survfit")){
    n_groups = length(data$n)
  }
  if(!is.null(data) & is.null(se_KM) & n_groups == 2){
    return(get_se_KM_from_data(data, t))
  }
  else if(length(se_KM) == 1 & (n_groups == 1| n_groups == 0) & is.numeric(se_KM)){
    if(n_groups == 0){warning("No group variable was detected in the data frame. Assuming that there is only one group.")}
    return(c(get_se_KM_from_data(data,t), se_KM))
  } else if(length(se_KM) == 2 & is.null(data) & is.numeric(se_KM)){
    if(!is.null(t)){
      warning("If you want to compare two numeric values of survival
            the time variable t is not used.")
    }
    return(se_KM)
  } else{
    stop("Invalid input. There are three possibilities for data input,
        (1) a data frame or a survfit object 'data' with two groups and null variables 'surv_KM' and
        'se_KM', (2) a data frame or a survfit object 'data' with one group and numeric vectors 'surv_KM' and
        'se_KM' of length 1, (3) a null variable 'data' and numeric vectors 'surv_KM'
        and 'se_KM' of length 2.")
  }
}


#' Square root of Greenwood's formula
#'
#' This function returns the square root of the quotient of the Kaplan-Meier
#' estimate of the survival rate by its standard error. This quotient can be
#' calculated by a formula known as Greenwood's formula. This helper function is
#' simply introduced for better readability when comparing the implemented
#' functions to the ones from Klein et al.
#'
#' @param surv_KM a numeric vector containing survival rates for two groups
#' @param se_KM a numeric vector containing standard errors for two groups
#'
#' @return A numeric vector containing the square root of Greenwood's formula for
#'  the different groups.
#'
#' @examples
#' get_sigma_KM(c(0.8, 0.9), c(0.11, 0.09))
#'
#' @export
get_sigma_KM <- function(surv_KM, se_KM){
  return(se_KM/surv_KM)
}

#' Function for generating a survfit object.
#'
#' This function returns a \code{\link[survival]{survfit}} object that has been generated from a
#' survival data set.
#'
#' @param data a data frame that has the variables \code{time} specifying time to
#' event or time to censoring, \code{status} specifying the censoring status and
#' \code{group} specifying the different groups to be compared.
#'
#' @return A \code{\link[survival]{survfit}} object.
#'
#' @examples
#' data(exp_surv)
#' get_survfit(exp_surv)
#'
#' @export
get_survfit <- function(data, group = NULL, time=NULL, status=NULL){
  # Better idea: Set time = time, status = status, group= group as a standard, then everything gets easier :-)
  # If data is of type survfit or null, don't change anything.
  if(methods::is(data, "survfit") | is.null(data)){
    return(data)
  }
  # Check whether variables time and status have been specified
  time_exists = "time" %in% names(data)
  status_exists = "status" %in% names(data)
  if(!(time_exists && status_exists)) {
    stop("Data set must either contain the variables time and status,
         or the corresponding variables must be specified via the time= and group= options.")
  }
  # check for numeric time variable.
  if (!is.numeric(data$time)) {
    stop("Variable time must be a numeric vector.")
  }
  # Split by group, in case such a variable exists.
  if("group" %in% names(data)){
    return(survival::survfit(formula = survival::Surv(time, status) ~ group, data=data))
  }
  else if(!is.null(group)){
    return(survival::survfit(formula = survival::Surv(time, status) ~ group, data=data))
  } else{
    return(survival::survfit(formula = survival::Surv(time, status) ~ group, data=data))
  }
}


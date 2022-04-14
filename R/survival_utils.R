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
get_surv_KM_from_data <- function(data, t, ...){
  # TODO: Ist das überhaupt in Ordnung, dass ich hier extend = TRUE gewählt habe?
  if(!methods::is(data, "survfit")){
    data_survfit <- get_survfit(data, ...)
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
get_surv_KM <- function(data = NULL, surv_KM = NULL, t = NULL, ...){
  n_groups = 0
  if(methods::is(data, "data.frame")){
    n_groups = length(levels(as.factor(data$group)))
  }
  else if(methods::is(data, "survfit")){
    n_groups = length(data$n)
  }
  if(!is.null(data) & is.null(surv_KM) & n_groups == 2){
    return(get_surv_KM_from_data(data, t, ...))
  }
  else if(length(surv_KM) == 1 & (n_groups == 1 | n_groups == 0) & is.numeric(surv_KM)){
    if(n_groups == 0){warning("No group variable was detected in the data frame. Assuming that there is only one group.")}
    return(c(get_surv_KM_from_data(data,t), surv_KM, ...))
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
get_se_KM <- function(data = NULL, se_KM = NULL, t = NULL, ...){
  n_groups = 0
  if(methods::is(data, "data.frame")){
    n_groups = length(levels(as.factor(data$group)))
  }
  else if(methods::is(data, "survfit")){
    n_groups = length(data$n)
  }
  if(!is.null(data) & is.null(se_KM) & n_groups == 2){
    return(get_se_KM_from_data(data, t, ...))
  }
  else if(length(se_KM) == 1 & (n_groups == 1| n_groups == 0) & is.numeric(se_KM)){
    if(n_groups == 0){warning("No group variable was detected in the data frame. Assuming that there is only one group.")}
    return(c(get_se_KM_from_data(data, t, ...), se_KM))
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

get_survfit <- function(data, time = time, status = status, group = group){
  time = rlang::enquo(time)
  if (rlang::quo_is_null(time)) {
    warning("You supplied time=NULL as an argument. Assuming you meant the
            default value time=time.")
    time = rlang::sym("time")
  } else if (rlang::quo_is_symbol(time)) {
    time = rlang::get_expr(time)
  } else {
    stop(paste("Expected symbol but found time=", class(rlang::get_expr(time))))
  }
  status = rlang::enquo(status)
  if (rlang::quo_is_null(status)) {
    warning("You supplied status=NULL as an argument. Assuming you meant the
            default value status=status.")
    status = rlang::sym("status")
  } else if (rlang::quo_is_symbol(status)) {
    status = rlang::get_expr(status)
  } else {
    stop(paste("Expected symbol but found status=", class(rlang::get_expr(status))))
  }
  group = rlang::enquo(group)
  if (rlang::quo_is_null(group)) {
    warning("You supplied group=NULL as an argument. Assuming you meant the
            default value group=1, i.e. no grouping of survival data.")
    group = 1
  } else if (rlang::quo_is_symbol(group)) {
    group = rlang::get_expr(group)
  } else if (rlang::get_expr(group) == 1){
    group = 1
  } else {
    stop(paste("Expected symbol or 1 or NULL but found group=", class(rlang::get_expr(status))))
  }

  # If data is of type survfit or null, don't change anything.

  if(methods::is(data, "survfit") | is.null(data)){
    return(data)
  }

  # Check whether data has variables time and status
  time_name = rlang::inject(deparse(substitute(!!time)))
  print(paste("Time name =",time_name))
  if(!(as.character(time_name) %in% names(data))) {
    stop("The time variable ", time_name, " is not in your data set.")
  }
  # Check whether data has variables time and status
  status_name = rlang::inject(deparse(substitute(!!status)))
  if(!(as.character(status_name) %in% names(data))) {
    stop("The status variable ", status_name, " is not in your data set.")
  }
  # Split by group, in case such a variable exists.
  if(as.character(group) %in% names(data) | group == 1){
    return(
    rlang::inject(
      survival::survfit(formula = survival::Surv(!!time, !!status) ~ !!group, data)
    )
    )
  }
  else{
    stop("The group variable ", group, " is not in your data set.")
  }
}


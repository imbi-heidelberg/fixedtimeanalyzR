# ------------------------------------------------------------------------------
# STATISTICAL TESTS FOR COMPARING SURVIVAL CURVES


#' Naive test for difference of survival rates
#'
#' This function tests two survival rates at a fixed point in time for equality.
#' This test can also be used in order to compare the survival rate of a data
#'  frame to some given numerical value. Therefore, there are three
#' possibilities for data input, (1) \code{data} with two groups and
#' null variables \code{surv_KM} and \code{se_KM}, (2) \code{data} with one group and two
#' numeric vectors \code{surv_KM} and \code{se_KM} of length 1, (3) a null variable \code{data} and
#' two numeric vectors \code{surv_KM} and \code{se_KM} of length 2. The test uses a naive
#'  test statistic for comparing the survival rates of two survival curves at a
#'  fixed point in time. This statistic is an implementation of the
#' test statistic \eqn{X_1^2}{X1^2} in Klein et al.'s paper [1]. The test is a
#' chi-squared test with one degree of freedom.
#'
#' @references [1]  Klein, J.P., Logan, B., Harhoff, M. and Andersen, P.K. (2007),
#' \href{https://doi.org/10.1002/sim.2864}{Analyzing survival curves at a fixed point in time}.
#' Statist. Med., 26: 4505-4519.
#'
#' @param data either a null object or a \code{\link[survival]{survfit}} object or a data frame that has
#' the variables \code{time} specifying time to event or time to censoring,
#' \code{status} specifying the censoring status and \code{group} specifying the
#'  different groups to be compared.
#' @param surv_KM either a null object or a numeric vector containing survival
#' rates for one or two groups
#' @param se_KM either a null object or a numeric vector containing standard
#' errors for one or two groups
#' @param t the point in time at which the survival rate will be calculated
#'
#' @return A list of class \code{\link[EnvStats:htest.object]{htest}} containing information about the
#' p-value, the value of the test statistic, survival rates etc.
#'
#' @examples
#' data(exp_surv)
#' naive.test(exp_surv, t=1)
#'
#' library(survival)
#' exp_fit <- survfit(formula = Surv(time, status) ~ group, data=exp_surv)
#' naive.test(exp_fit, t=1)
#'
#' library(dplyr)
#' exp_surv %>% filter(group == 1) -> exp_surv_grp1
#' naive.test(exp_surv_grp1, surv_KM=0.5, se_KM=0.08, t=1)
#'
#' naive.test(surv_KM=c(0.5, 0.6), se_KM=c(0.1,0.1))
#'
#' @export
naive.test <- function(data = NULL, surv_KM = NULL, se_KM = NULL, t = NULL,
                      time= time, status= status, group = group) {
  time = rlang::enquo(time)
  status = rlang::enquo(status)
  group = rlang::enquo(group)
  # Save data name
  res <- list()
  res$data.name   <- sprintf(deparse(substitute(data)))
  # Convert data to survfit object - nothing happens when data is already survfit
  data = rlang::inject(get_survfit(data, !!time, !!status, !!group))
  # Calculation and/or concatenation of survival data
  surv_KM = get_surv_KM(data, surv_KM, t)
  se_KM = get_se_KM(data, se_KM, t)


  # Construction of the htest results object
  class(res) <- "htest"
  res$null.value  <- c(0)
  names(res$null.value) <- c("absolute survival rate difference |group 1 - group 2|")
  res$alternative <- "greater" # A chi-squared test is always a one-sided test
  res$method <- "Naive test for equality of survival rates at a fixed point in time.

  This test is using the naive test statistic X1^2 from Klein et al.'s paper
  'Analyzing survival curves at a fixed point in time', published in Stat. Med.,
  26 (2007). It is a chi-squared test with one degree of freedom."
  if(!is.null(t)){
    res$parameters  <- c(t)
    names(res$parameters) <- c("time")
  }

  # Survival rates and their difference
  if(!is.null(data$strata)){
    grp_names <- names(data$strata)
  } else if(!is.null(data)){
    grp_names <- c("data", "numeric input")
  } else{
    grp_names <- c("numeric input 1", "numeric input 2")
  }
  diff_surv_KM <- abs(surv_KM[1] - surv_KM[2])
  res$estimate <- c(surv_KM[1], se_KM[1], surv_KM[2], se_KM[2], diff_surv_KM)
  names(res$estimate) <- c(paste("Survival rate of ", grp_names[1]),
                           paste("Standard error of ", grp_names[1]),
                           paste("Survival rate of ", grp_names[2]),
                           paste("Standard error of ", grp_names[2]),
                           "Absolute difference of rates")

  # Test statistic
  statistic <- naive.t(surv_KM, se_KM)
  res$statistic <- statistic
  names(res$statistic) <- "naive statistic X1^2"

  # p-value, probability of naive.t < t
  res$p.value <- 1 - stats::pchisq(statistic, df = 1)
  return(res)
}

#' Log-transformed test for difference of survival rates
#'
#' This function tests two survival rates at a fixed point in time for equality.
#' This test can also be used in order to compare the survival rate of a data
#'  frame to some given numerical value. Therefore, there are three
#' possibilities for data input, (1) \code{data} with two groups and
#' null variables \code{surv_KM} and \code{se_KM}, (2) \code{data} with one group and two
#' numeric vectors \code{surv_KM} and \code{se_KM} of length 1, (3) a null variable \code{data} and
#' two numeric vectors \code{surv_KM} and \code{se_KM} of length 2. The test uses a log-transformed
#'  test statistic for comparing the survival rates of two survival curves at a
#'  fixed point in time. This statistic is an implementation of the
#' test statistic \eqn{X_2^2}{X2^2} in Klein et al.'s paper [1]. The test is a
#' chi-squared test with one degree of freedom.
#'
#' @references [1]  Klein, J.P., Logan, B., Harhoff, M. and Andersen, P.K. (2007),
#' \href{https://doi.org/10.1002/sim.2864}{Analyzing survival curves at a fixed point in time}.
#' Statist. Med., 26: 4505-4519.
#'
#' @param data either a null object or a \code{\link[survival]{survfit}} object or a data frame that has
#' the variables \code{time} specifying time to event or time to censoring,
#' \code{status} specifying the censoring status and \code{group} specifying the
#'  different groups to be compared.
#' @param surv_KM either a null object or a numeric vector containing survival
#' rates for one or two groups
#' @param se_KM either a null object or a numeric vector containing standard
#' errors for one or two groups
#' @param t the point in time at which the survival rate will be calculated
#'
#' @return A list of class \code{\link[EnvStats:htest.object]{htest}} containing information about the
#' p-value, the value of the test statistic, survival rates etc.
#'
#' @examples
#' data(exp_surv)
#' logtra.test(exp_surv, t=1)
#'
#' library(survival)
#' exp_fit <- survfit(formula = Surv(time, status) ~ group, data=exp_surv)
#' logtra.test(exp_fit, t=1)
#'
#' library(dplyr)
#' exp_surv %>% filter(group == 1) -> exp_surv_grp1
#' logtra.test(exp_surv_grp1, surv_KM=0.5, se_KM=0.08, t=1)
#'
#' logtra.test(surv_KM=c(0.5, 0.6), se_KM=c(0.1,0.1))
#'
#' @export
logtra.test <- function(data = NULL, surv_KM = NULL, se_KM = NULL, t = NULL,
                        time= time, status= status, group = group) {
  time = rlang::enquo(time)
  status = rlang::enquo(status)
  group = rlang::enquo(group)
  # Save data name
  res <- list()
  res$data.name   <- sprintf(deparse(substitute(data)))
  # Convert data to survfit object - nothing happens when data is already survfit
  data = rlang::inject(get_survfit(data, !!time, !!status, !!group))

  # Calculation and/or concatenation of survival data
  surv_KM = get_surv_KM(data, surv_KM, t)
  se_KM = get_se_KM(data, se_KM, t)

  # Construction of the htest results object

  class(res) <- "htest"
  res$null.value  <- c(0)
  names(res$null.value) <- c("absolute survival rate difference |group 1 - group 2|")
  res$alternative <- "greater" # A chi-squared test is always a one-sided test
  res$method <- "Log-transformed test for equality of survival rates at a fixed point in time.

  This test is using the log-transformed test statistic X2^2 from Klein et al.'s paper
  'Analyzing survival curves at a fixed point in time', published in Stat. Med.,
  26 (2007). It is a chi-squared test with one degree of freedom."

  if(!is.null(t)){
    res$parameters  <- c(t)
    names(res$parameters) <- c("time")
  }

  # Survival rates and their difference
  grp_names <- names(data$strata)
  diff_surv_KM <- abs(surv_KM[1] - surv_KM[2])
  res$estimate <- c(surv_KM[1], se_KM[1], surv_KM[2], se_KM[2], diff_surv_KM)
  names(res$estimate) <- c("Survival rate of group 1",
                           "Standard error of group 1",
                           "Survival rate of group 2",
                           "Standard error of group 2",
                           "Absolute difference of rates")

  # Test statistic
  statistic <- logtra.t(surv_KM, se_KM)
  res$statistic <- statistic
  names(res$statistic) <- "log-transformed statistic X2^2"

  # p-value, probability of logtra.t < t
  res$p.value <- 1 - stats::pchisq(statistic, df = 1)
  return(res)
}

#' Clog-transformed test for difference of survival rates
#'
#' This function tests two survival rates at a fixed point in time for equality.
#' This test can also be used in order to compare the survival rate of a data
#'  frame to some given numerical value. Therefore, there are three
#' possibilities for data input, (1) \code{data} with two groups and
#' null variables \code{surv_KM} and \code{se_KM}, (2) \code{data} with one group and two
#' numeric vectors \code{surv_KM} and \code{se_KM} of length 1, (3) a null variable \code{data} and
#' two numeric vectors \code{surv_KM} and \code{se_KM} of length 2. The test uses a clog-transformed
#'  test statistic for comparing the survival rates of two survival curves at a
#'  fixed point in time. This statistic is an implementation of the
#' test statistic \eqn{X_3^2}{X3^2} in Klein et al.'s paper [1]. The test is a
#' chi-squared test with one degree of freedom.
#'
#' @references [1]  Klein, J.P., Logan, B., Harhoff, M. and Andersen, P.K. (2007),
#' \href{https://doi.org/10.1002/sim.2864}{Analyzing survival curves at a fixed point in time}.
#' Statist. Med., 26: 4505-4519.
#'
#' @param data either a null object or a \code{\link[survival]{survfit}} object or a data frame that has
#' the variables \code{time} specifying time to event or time to censoring,
#' \code{status} specifying the censoring status and \code{group} specifying the
#'  different groups to be compared.
#' @param surv_KM either a null object or a numeric vector containing survival
#' rates for one or two groups
#' @param se_KM either a null object or a numeric vector containing standard
#' errors for one or two groups
#' @param t the point in time at which the survival rate will be calculated
#'
#' @return A list of class \code{\link[EnvStats:htest.object]{htest}} containing information about the
#' p-value, the value of the test statistic, survival rates etc.
#'
#' @examples
#' data(exp_surv)
#' clog.test(exp_surv, t=1)
#'
#' library(survival)
#' exp_fit <- survfit(formula = Surv(time, status) ~ group, data=exp_surv)
#' clog.test(exp_fit, t=1)
#'
#' library(dplyr)
#' exp_surv %>% filter(group == 1) -> exp_surv_grp1
#' clog.test(exp_surv_grp1, surv_KM=0.5, se_KM=0.08, t=1)
#'
#' clog.test(surv_KM=c(0.5, 0.6), se_KM=c(0.1,0.1))
#'
#' @export
clog.test <- function(data = NULL, surv_KM = NULL, se_KM = NULL, t = NULL,
                      time= time, status= status, group = group) {
  time = rlang::enquo(time)
  status = rlang::enquo(status)
  group = rlang::enquo(group)
  # Save data name
  res <- list()
  res$data.name   <- sprintf(deparse(substitute(data)))
  # Convert data to survfit object - nothing happens when data is already survfit
  data = rlang::inject(get_survfit(data, !!time, !!status, !!group))
  # Calculation and/or concatenation of survival data
  surv_KM = get_surv_KM(data, surv_KM, t)
  se_KM = get_se_KM(data, se_KM, t)


  # Construction of the htest results object
  class(res) <- "htest"
  res$null.value  <- c(0)
  names(res$null.value) <- c("absolute survival rate difference |group 1 - group 2|")
  res$alternative <- "greater" # A chi-squared test is always a one-sided test
  res$method <- "Clog-transformed test for equality of survival rates at a fixed point in time.

  This test is using the clog-transformed test statistic X2^2 from Klein et al.'s paper
  'Analyzing survival curves at a fixed point in time', published in Stat. Med.,
  26 (2007). It is a chi-squared test with one degree of freedom."
  if(!is.null(t)){
    res$parameters  <- c(t)
    names(res$parameters) <- c("time")
  }

  # Survival rates and their difference
  grp_names <- names(data$strata)
  diff_surv_KM <- abs(surv_KM[1] - surv_KM[2])
  res$estimate <- c(surv_KM[1], se_KM[1], surv_KM[2], se_KM[2], diff_surv_KM)
  names(res$estimate) <- c(paste("Survival rate of", grp_names[1]),
                           paste("Standard error of", grp_names[1]),
                           paste("Survival rate of", grp_names[2]),
                           paste("Standard error of", grp_names[2]),
                           "Absolute difference of rates")

  # Test statistic
  statistic <- clog.t(surv_KM, se_KM)
  res$statistic <- statistic
  names(res$statistic) <- "clog-transformed statistic X3^2"

  # p-value, probability of clog.t < t
  res$p.value <- 1 - stats::pchisq(statistic, df = 1)
  return(res)
}

#' Arc sine and square root-transformed test for difference of survival rates
#'
#' This function tests two survival rates at a fixed point in time for equality.
#' This test can also be used in order to compare the survival rate of a data
#'  frame to some given numerical value. Therefore, there are three
#' possibilities for data input, (1) \code{data} with two groups and
#' null variables \code{surv_KM} and \code{se_KM}, (2) \code{data} with one group and two
#' numeric vectors \code{surv_KM} and \code{se_KM} of length 1, (3) a null variable \code{data} and
#' two numeric vectors \code{surv_KM} and \code{se_KM} of length 2. The test uses a test statistic
#' with an arc-sine square-root transformation for comparing the survival rates
#' of two survival curves at a fixed point in time. This statistic is an
#' implementation of the test statistic \eqn{X_4^2}{X4^2} in Klein et al.'s paper
#'  [1]. The test is a chi-squared test with one degree of freedom.
#'
#' @references [1]  Klein, J.P., Logan, B., Harhoff, M. and Andersen, P.K. (2007),
#' \href{https://doi.org/10.1002/sim.2864}{Analyzing survival curves at a fixed point in time}.
#' Statist. Med., 26: 4505-4519.
#'
#' @param data either a null object or a \code{\link[survival]{survfit}} object or a data frame that has
#' the variables \code{time} specifying time to event or time to censoring,
#' \code{status} specifying the censoring status and \code{group} specifying the
#'  different groups to be compared.
#' @param surv_KM either a null object or a numeric vector containing survival
#' rates for one or two groups
#' @param se_KM either a null object or a numeric vector containing standard
#' errors for one or two groups
#' @param t the point in time at which the survival rate will be calculated
#'
#' @return A list of class \code{\link[EnvStats:htest.object]{htest}} containing information about the
#' p-value, the value of the test statistic, survival rates etc.
#'
#' @examples
#' data(exp_surv)
#' asinsqrt.test(exp_surv, t=1)
#'
#' library(survival)
#' exp_fit <- survfit(formula = Surv(time, status) ~ group, data=exp_surv)
#' asinsqrt.test(exp_fit, t=1)
#'
#' library(dplyr)
#' exp_surv %>% filter(group == 1) -> exp_surv_grp1
#' asinsqrt.test(exp_surv_grp1, surv_KM=0.5, se_KM=0.08, t=1)
#'
#' asinsqrt.test(surv_KM=c(0.5, 0.6), se_KM=c(0.1,0.1))
#'
#' @export
asinsqrt.test <- function(data = NULL, surv_KM = NULL, se_KM = NULL, t = NULL,
                          time= time, status= status, group = group) {
  time = rlang::enquo(time)
  status = rlang::enquo(status)
  group = rlang::enquo(group)
  # Save data name
  res <- list()
  res$data.name   <- sprintf(deparse(substitute(data)))
  # Convert data to survfit object - nothing happens when data is already survfit
  data = rlang::inject(get_survfit(data, !!time, !!status, !!group))
  # Calculation and/or concatenation of survival data
  surv_KM = get_surv_KM(data, surv_KM, t)
  se_KM = get_se_KM(data, se_KM, t)


  # Construction of the htest results object
  class(res) <- "htest"
  res$null.value  <- c(0)
  names(res$null.value) <- c("absolute survival rate difference |group 1 - group 2|")
  res$alternative <- "greater" # A chi-squared test is always a one-sided test
  res$method <- "Arc-sine and square-root transformed test for equality of survival rates at a fixed point in time.

  This test is test statistic X4^2 transformed using arc-sine and square-root
  from Klein et al.'s paper 'Analyzing survival curves at a fixed point in time',
  published in Stat. Med., 26 (2007). It is a chi-squared test with one degree of freedom."
  if(!is.null(t)){
    res$parameters  <- c(t)
    names(res$parameters) <- c("time")
  }

  # Survival rates and their difference
  diff_surv_KM <- abs(surv_KM[1] - surv_KM[2])
  res$estimate <- c(surv_KM[1], se_KM[1], surv_KM[2], se_KM[2], diff_surv_KM)
  names(res$estimate) <- c("Survival rate of group 1",
                           "Standard error of group 1",
                           "Survival rate of group 2",
                           "Standard error of group 2",
                           "Absolute difference of rates")

  # Test statistic
  statistic <- asinsqrt.t(surv_KM, se_KM)
  res$statistic <- statistic
  names(res$statistic) <- "asinsqrt-transformed statistic X4^2"

  # p-value, probability of asinsqrt.t < t
  res$p.value <- 1 - stats::pchisq(statistic, df = 1)
  return(res)
}

#' Logit-transformed test for difference of survival rates
#'
#' This function tests two survival rates at a fixed point in time for equality.
#' This test can also be used in order to compare the survival rate of a data
#'  frame to some given numerical value. Therefore, there are three
#' possibilities for data input, (1) \code{data} with two groups and
#' null variables \code{surv_KM} and \code{se_KM}, (2) \code{data} with one group and two
#' numeric vectors \code{surv_KM} and \code{se_KM} of length 1, (3) a null variable \code{data} and
#' two numeric vectors \code{surv_KM} and \code{se_KM} of length 2. The test uses a logit-transformed
#'  test statistic for comparing the survival rates of two survival curves at a
#'  fixed point in time. This statistic is an implementation of the
#' test statistic \eqn{X_5^2}{X5^2} in Klein et al.'s paper [1]. The test is a
#' chi-squared test with one degree of freedom.
#'
#' @references [1]  Klein, J.P., Logan, B., Harhoff, M. and Andersen, P.K. (2007),
#' \href{https://doi.org/10.1002/sim.2864}{Analyzing survival curves at a fixed point in time}.
#' Statist. Med., 26: 4505-4519.
#'
#' @param data either a null object or a \code{\link[survival]{survfit}} object or a data frame that has
#' the variables \code{time} specifying time to event or time to censoring,
#' \code{status} specifying the censoring status and \code{group} specifying the
#'  different groups to be compared.
#' @param surv_KM either a null object or a numeric vector containing survival
#' rates for one or two groups
#' @param se_KM either a null object or a numeric vector containing standard
#' errors for one or two groups
#' @param t the point in time at which the survival rate will be calculated
#'
#' @return A list of class \code{\link[EnvStats:htest.object]{htest}} containing information about the
#' p-value, the value of the test statistic, survival rates etc.
#'
#' @examples
#' data(exp_surv)
#' logit.test(exp_surv, t=1)
#'
#' library(survival)
#' exp_fit <- survfit(formula = Surv(time, status) ~ group, data=exp_surv)
#' logit.test(exp_fit, t=1)
#'
#' library(dplyr)
#' exp_surv %>% filter(group == 1) -> exp_surv_grp1
#' logit.test(exp_surv_grp1, surv_KM=0.5, se_KM=0.08, t=1)
#'
#' logit.test(surv_KM=c(0.5, 0.6), se_KM=c(0.1,0.1))
#'
#' @export
logit.test <- function(data = NULL, surv_KM = NULL, se_KM = NULL, t = NULL,
                       time= time, status= status, group = group) {
  time = rlang::enquo(time)
  status = rlang::enquo(status)
  group = rlang::enquo(group)
  # Save data name
  res <- list()
  res$data.name   <- sprintf(deparse(substitute(data)))
  # Convert data to survfit object - nothing happens when data is already survfit
  data = rlang::inject(get_survfit(data, !!time, !!status, !!group))
  # Calculation and/or concatenation of survival data
  surv_KM = get_surv_KM(data, surv_KM, t)
  se_KM = get_se_KM(data, se_KM, t)


  # Construction of the htest results object
  class(res) <- "htest"
  res$null.value  <- c(0)
  names(res$null.value) <- c("absolute survival rate difference |group 1 - group 2|")
  res$alternative <- "greater" # A chi-squared test is always a one-sided test
  res$method <- "logit-transformed test for equality of survival rates at a fixed point in time.

  This test is using the logit-transformed test statistic X5^2 from Klein et al.'s paper
  'Analyzing survival curves at a fixed point in time', published in Stat. Med.,
  26 (2007). It is a chi-squared test with one degree of freedom."
  if(!is.null(t)){
    res$parameters  <- c(t)
    names(res$parameters) <- c("time")
  }

  # Survival rates and their difference
  diff_surv_KM <- abs(surv_KM[1] - surv_KM[2])
  res$estimate <- c(surv_KM[1], se_KM[1], surv_KM[2], se_KM[2], diff_surv_KM)
  names(res$estimate) <- c("Survival rate of group 1",
                           "Standard error of group 1",
                           "Survival rate of group 2",
                           "Standard error of group 2",
                           "Absolute difference of rates")

  # Test statistic
  statistic <- logit.t(surv_KM, se_KM)
  res$statistic <- statistic
  names(res$statistic) <- "logit-transformed statistic X5^2"

  # p-value, probability of logit.t < t
  res$p.value <- 1 - stats::pchisq(statistic, df = 1)
  return(res)
}

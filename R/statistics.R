# ------------------------------------------------------------------------------
# TEST STATISTICS FOR COMPARING SURVIVAL CURVES


#' Naive test statistic.
#'
#' This function is a naive test statistic for comparing the survival rates of
#' two survival curves at a fixed point in time. It is an implementation of the
#' test statistic \eqn{X_1^2}{X1^2} in Klein et al.'s paper [1].
#'
#' @references [1]  Klein, J.P., Logan, B., Harhoff, M. and Andersen, P.K. (2007),
#' \href{https://doi.org/10.1002/sim.2864}{Analyzing survival curves at a fixed point in time}.
#' Statist. Med., 26: 4505-4519.
#'
#' @param surv_KM a numeric vector containing survival rates for two groups
#' @param se_KM a numeric vector containing standard errors for two groups
#'
#' @return A numeric scalar containing the value of the test statistic.
#'
#' @examples
#' data(exp_surv)
#'
#' surv_KM = get_surv_KM(exp_surv, t=1)
#' se_KM = get_se_KM(exp_surv, t=1)
#' naive.t(surv_KM, se_KM)
#'
#' @export
naive.t <- function(surv_KM, se_KM){
  return( (surv_KM[1]-surv_KM[2])**2/ (se_KM[1]**2 + se_KM[2]**2) )
}

#' Log-transformed test statistic.
#'
#' This function is a test statistic for comparing the survival rates of
#' two survival curves at a fixed point in time, using a logarithmic
#' transformation of the survival function. It is an implementation of the
#' test statistic \eqn{X_2^2}{X2^2} in Klein et al.'s paper [1].
#'
#' [1] - Klein, J.P., Logan, B., Harhoff, M. and Andersen, P.K. (2007),
#' \href{https://doi.org/10.1002/sim.2864}{Analyzing survival curves at a fixed point in time}.
#' Statist. Med., 26: 4505-4519.
#'
#' @param surv_KM a numeric vector containing survival rates for two groups
#' @param se_KM a numeric vector containing standard errors for two groups
#'
#' @return A numeric scalar containing the value of the test statistic.
#'
#' @examples
#' data(exp_surv)
#'
#' surv_KM = get_surv_KM(exp_surv, t=1)
#' se_KM = get_se_KM(exp_surv, t=1)
#' logtra.t(surv_KM, se_KM)
#'
#' @export
logtra.t <- function(surv_KM, se_KM){
  sigma_KM = get_sigma_KM(surv_KM, se_KM)
  return( (log(surv_KM[1])-log(surv_KM[2]))**2/ (sigma_KM[1]**2 + sigma_KM[2]**2) )
}

#' Clog-transformed test statistic.
#'
#' This function is a test statistic for comparing the survival rates of
#' two survival curves at a fixed point in time, using a log(-log)-
#' transformation of the survival function. It is an implementation of the
#' test statistic \eqn{X_3^2}{X3^2} in Klein et al.'s paper [1].
#'
#' @references [1]  Klein, J.P., Logan, B., Harhoff, M. and Andersen, P.K. (2007),
#' \href{https://doi.org/10.1002/sim.2864}{Analyzing survival curves at a fixed point in time}.
#' Statist. Med., 26: 4505-4519.
#'
#' @param surv_KM a numeric vector containing survival rates for two groups
#' @param se_KM a numeric vector containing standard errors for two groups
#'
#' @return A numeric scalar containing the value of the test statistic.
#'
#' @examples
#' data(exp_surv)
#'
#' surv_KM = get_surv_KM(exp_surv, t=1)
#' se_KM = get_se_KM(exp_surv, t=1)
#' clog.t(surv_KM, se_KM)
#'
#' @export
clog.t <- function(surv_KM, se_KM){
  sigma_KM = get_sigma_KM(surv_KM, se_KM)
  # TODO: Warning! Due to the fact that log(1)=0, this statistic can return NaN, in case that surv_KM equals 1.
  numerator = (log(-log(surv_KM[1]))-log(-log(surv_KM[2])))**2
  denominator = ( sigma_KM[1]/log(surv_KM[1]) )**2+( sigma_KM[2]/log(surv_KM[2]) )**2
  return(numerator/denominator)
}

#' Arc sine and square root-transformed test statistic.
#'
#' This function is a test statistic for comparing the survival rates of
#' two survival curves at a fixed point in time, using an arc-sine square-root
#' transformation of the survival function. It is an implementation of the
#' test statistic \eqn{X_4^2}{X4^2} in Klein et al.'s paper [1].
#'
#' @references [1]  Klein, J.P., Logan, B., Harhoff, M. and Andersen, P.K. (2007),
#' \href{https://doi.org/10.1002/sim.2864}{Analyzing survival curves at a fixed point in time}.
#' Statist. Med., 26: 4505-4519.
#'
#' @param surv_KM a numeric vector containing survival rates for two groups
#' @param se_KM a numeric vector containing standard errors for two groups
#'
#' @return A numeric scalar containing the value of the test statistic.
#'
#' @examples
#' data(exp_surv)
#'
#' surv_KM = get_surv_KM(exp_surv, t=1)
#' se_KM = get_se_KM(exp_surv, t=1)
#' asinsqrt.t(surv_KM, se_KM)
#'
#' @export
asinsqrt.t <- function(surv_KM, se_KM){
  sigma_KM = get_sigma_KM(surv_KM, se_KM)
  nu = c(0.0,0.0)
  for(k in (1:2)){
    nu[k] = surv_KM[k]*sigma_KM[k]**2/( 4*(1-surv_KM[k]) )
  }
  return( ( asin(sqrt(surv_KM[1]))-asin(sqrt(surv_KM[2])) )**2/ (nu[1]+nu[2]) )
}

#' Logit-transformed test statistic.
#'
#' This function is a test statistic for comparing the survival rates of
#' two survival curves at a fixed point in time, using a logit-
#' transformation of the survival function. It is an implementation of the
#' test statistic \eqn{X_5^2}{X5^2} in Klein et al.'s paper [1].
#'
#' @references [1]  Klein, J.P., Logan, B., Harhoff, M. and Andersen, P.K. (2007),
#' \href{https://doi.org/10.1002/sim.2864}{Analyzing survival curves at a fixed point in time}.
#' Statist. Med., 26: 4505-4519.
#'
#' @param surv_KM a numeric vector containing survival rates for two groups
#' @param se_KM a numeric vector containing standard errors for two groups
#'
#' @return A numeric scalar containing the value of the test statistic.
#'
#' @examples
#' data(exp_surv)
#'
#' surv_KM = get_surv_KM(exp_surv, t=1)
#' se_KM = get_se_KM(exp_surv, t=1)
#' logit.t(surv_KM, se_KM)
#'
#' @export
logit.t <- function(surv_KM, se_KM){
  sigma_KM = get_sigma_KM(surv_KM, se_KM)
  numerator = ( log(surv_KM[1]/(1-surv_KM[1])) - log(surv_KM[2]/(1-surv_KM[2])) )**2
  denominator = sigma_KM[1]**2/(1-surv_KM[1])**2 + sigma_KM[2]**2/(1-surv_KM[2])**2
  return(numerator/denominator)
}


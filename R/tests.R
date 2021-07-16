library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(tidyverse)
library(magrittr)








X1.test <- function(surv_KM, se_KM, conf.level = 0.95){
  statistic <- X1(surv_KM, se_KM)
  result <- statistic > qchisq(conf.level, 1)
  l = c(statistic, result, "Naive test X1")
  names(l) <- c("statistic", "result", "method")

  return(l)
}

#-------------------------------------------------------------------------------
# log-trafo statistic
X2 <- function(surv_KM, se_KM){
  sigma_KM = get_sigma_KM(surv_KM, se_KM)
  return( (log(surv_KM[1])-log(surv_KM[2]))**2/ (sigma_KM[1]**2 + sigma_KM[2]**2) )
}

X2.test <- function(surv_KM, se_KM, conf.level = 0.95){
  statistic <- X2(surv_KM, se_KM)
  result <- statistic > qchisq(conf.level, 1)
  l = c(statistic, result, "log-transformed statistic X2")
  names(l) <- c("statistic", "result", "method")

  return(l)
}

#-------------------------------------------------------------------------------
# clog-trafo statistic
X3 <- function(surv_KM, se_KM){
  sigma_KM = get_sigma_KM(surv_KM, se_KM)
  # TODO: Warning! Due to the fact that log(1)=0, this statistic can return NaN, in case that surv_KM equals 1.
  numerator = (log(-log(surv_KM[1]))-log(-log(surv_KM[2])))**2
  denominator = ( sigma_KM[1]/log(surv_KM[1]) )**2+( sigma_KM[2]/log(surv_KM[2]) )**2
  return(numerator/denominator)
}

X3.test <- function(surv_KM, se_KM, conf.level = 0.95){
  statistic <- X3(surv_KM, se_KM)
  result <- statistic > qchisq(conf.level, 1)
  l = c(statistic, result, "clog-transformed statistic X3")
  names(l) <- c("statistic", "result", "method")

  return(l)
}

#-------------------------------------------------------------------------------
# arcsin-sqrt-trafo statistic
X4 <- function(surv_KM, se_KM){
  sigma_KM = get_sigma_KM(surv_KM, se_KM)
  nu = c(0.0,0.0)
  for(k in (1:2)){
    nu[k] = surv_KM[k]*sigma_KM[k]**2/( 4*(1-surv_KM[k]) )
  }
  return( ( asin(sqrt(surv_KM[1]))-asin(sqrt(surv_KM[2])) )**2/ (nu[1]+nu[2]) )
}

X4.test <- function(surv_KM, se_KM, conf.level = 0.95){
  statistic <- X4(surv_KM, se_KM)
  result <- statistic > qchisq(conf.level, 1)
  l = c(statistic, result, "arcsine-square-root-transformed statistic X4")
  names(l) <- c("statistic", "result", "method")

  return(l)
}

#-------------------------------------------------------------------------------
# logit-trafo statistic
X5 <- function(surv_KM, se_KM){
  sigma_KM = get_sigma_KM(surv_KM, se_KM)
  numerator = ( log(surv_KM[1]/(1-surv_KM[1])) - log(surv_KM[2]/(1-surv_KM[2])) )**2
  denominator = sigma_KM[1]**2/(1-surv_KM[1])**2 + sigma_KM[2]**2/(1-surv_KM[2])**2
  return(numerator/denominator)
}

X5.test <- function(surv_KM, se_KM, conf.level = 0.95){
  statistic <- X5(surv_KM, se_KM)
  result <- statistic > qchisq(conf.level, 1)
  l = c(statistic, result, "logit-transformed statistic X5")
  names(l) <- c("statistic", "result", "method")

  return(l)
}

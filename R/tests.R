library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(tidyverse)
library(magrittr)

#-------------------------------------------------------------------------------
# GENERATING RANDOM SAMPLE AND PLOTS FOR SURVIVAL ANALYSIS

# ------------------------------------------------------------------------------
# Generating a sample of two exponential distributions with exponential censoring times
#
# N -- is the number of all patients
# d -- is the factor between the first and second treatment group, i.e. n2 = d*n1, where N = n1 + n2
# or -- is the odds ratio between the two treatments groups
# p1 -- is the probability of having an event in treatment group one at time t=1.
generate_exp_sample <- function(N, d, or, p1){
  n1 = N/d # Number of patients in the first sample
  n2 = N - n1 # Number of patients in the second sample

  k1 = log(1/p1) # Structure parameter of the first exponential distribution
  k2 = log(or*(1-p1)/p1+1) # Structure parameter of the second exponential distribution

  # Generating the samples from the two data frames


  exp_sample = rbind(
    tibble(time=rexp(n1, k1), cens_time=rexp(n1, 0.1), treatment=rep(1, n1)),
    tibble(time=rexp(n2, k2), cens_time=rexp(n2, 0.1), treatment=rep(2, n2)))
  # status = TRUE if event has been observed, = FALSE if censored
  exp_sample %>% mutate(status = (cens_time > time)) %>%
    mutate(time = ifelse(status, time, cens_time)) %>%
    mutate(status = ifelse(status, 1, 0)) -> exp_sample

  return(exp_sample)
}

# ------------------------------------------------------------------------------
# Generating a plot from the sample, together with a exponential curve plot for indicating the distribution of
# the first treatment group.
#
# exp_sample -- is a sample with columns time, status and treatment
# p1 -- is the probability of having an event in treatment group 1 at time t = 1 (used as a shape parameter for
#       the exponential function)
generate_plot <- function(exp_sample, p1){
  # Computing and visualizing survival data
  exp_fit <- survfit(Surv(time, status) ~ treatment, data=exp_sample)
  plot = autoplot(exp_fit) + stat_function(fun=function(x){1-pexp(x, log(1/p1))}, geom="line", colour="orange")
  return(plot)
}

# What does this look like?
N = 90  # in c(30, 60, 90, 120, 150, 300)  N = n1 + n2 is the number of patients in both samples together
d = 2 #  in c(2, 3)  We either have n1=n2 or n2=2*n1, i.e. n1 = N/d
or = 2 # in c(1, 1.5, 2) Odds ratio between the two samples
p1 = 0.5 # in c(0.25,0.50,0.75) Probability of survival at time 1 in the first sample
sample = generate_exp_sample(N, d, or, p1)
plot = generate_plot(sample, p1)
plot


# ------------------------------------------------------------------------------
# DIFFERENT STATISTICS FOR COMPARING SURVIVAL CURVES

#-------------------------------------------------------------------------------
# Return the Kaplan-Meier estimate of the survival rate at a fixed time.
#
# sample -- a dataset with columns status, time and treatment
# t -- a number specifying the time at which survival rate shall be calculated
get_surv_KM <- function(sample, t){
  # TODO: Ist das überhaupt in Ordnung, dass ich hier extend = TRUE gewählt habe?
  surv_KM = summary(survfit(formula = Surv(time, status) ~ treatment, data=sample), times=t, extend = TRUE)$surv
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


# ------------------------------------------------------------------------------
# Naive statistic
X1 <- function(surv_KM, se_KM){
  return( (surv_KM[1]-surv_KM[2])**2/ (se_KM[1]**2 + se_KM[2]**2) )
}

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

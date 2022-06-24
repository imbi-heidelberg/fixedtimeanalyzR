# In this code test, we test the output of the statistical tests.
# As a reference, the two packages ComparisonSurv and bpcp are used.
# The test statistics...
# ... naive.test, logtra.test and clog.test are compared to results from bpcp
# ... clog.test, asinsqrt.test and logit.test are compared to results form ComparisonSurv.
#
#
# ComparisonSurv has two issues:
# (1) apparently, the statistics for naive and logtra were switched, so we omit these.
# (2) the p-values are only given up to the fourth digit.
# bpcp has one issue:
# (3) the statics for asinsqrt and logit were not implemented.

context("output of statistical tests")

library(dplyr)
library(ComparisonSurv)
library(bpcp)
library(survival)

# Preparations of data frames
data(exp_surv)
data(cancer, package = "survival")
exp_surv %>% mutate(group = group - 1) -> exp_surv # ComparisonSurv only accepts 0 and 1 as group variables
aml %>% rename(group = x) %>% mutate(group = factor(group, levels=c("Maintained","Nonmaintained"),labels=c(1,0)))-> amll # fixedtimeanalyzR only accepts group variable named group
amll_fit = survfit(formula = Surv(time, status) ~ group, data=amll)


# Function to apply all fixedtimeanaylzR test statistics to data.
apply_all_tests = function(data=NULL, surv_KM=NULL, se_KM=NULL, t=NULL) {
  naive = naive.test(data, surv_KM, se_KM, t = t)
  logtra = logtra.test(data, surv_KM, se_KM, t = t)
  clog = clog.test(data, surv_KM, se_KM, t = t)
  asinsqrt = asinsqrt.test(data, surv_KM, se_KM, t = t)
  logit = logit.test(data, surv_KM, se_KM, t = t)
  p.values = c(naive$p.value,
               logtra$p.value,
               clog$p.value,
               asinsqrt$p.value,
               logit$p.value)
  statistics = c(
    naive$statistic,
    logtra$statistic,
    clog$statistic,
    asinsqrt$statistic,
    logit$statistic
  )
  data.names = c (
    naive$data.name,
    logtra$data.name,
    clog$data.name,
    asinsqrt$data.name,
    logit$data.name
  )
  # Remove names, as references use different names
  names(p.values) = NULL
  names(statistics) = NULL
  return(list(
    p.values = p.values,
    statistics = statistics,
    data.names = data.names
  ))
}
# Function to apply all bpcp test statistics to data.
apply_bpcp_tests = function(data, t) {
  naive = bpcp::fixtdiff(
    data$time,
    data$status,
    data$group,
    trans = "identity",
    varpooled = FALSE,
    testtime = t
  )
  logtra = bpcp::fixtdiff(
    data$time,
    data$status,
    data$group,
    trans = "log",
    varpooled = FALSE,
    testtime = t
  )
  clog = bpcp::fixtdiff(
    data$time,
    data$status,
    data$group,
    trans = "cloglog",
    varpooled = FALSE,
    testtime = t
  )
  p.values = c(naive$p2, logtra$p2, clog$p2)
  names(p.values) = NULL
  return(list(p.values = p.values))

}
# Wrapper function to apply ComparisonSurv tests to data.
apply_ComparisonSurv_tests = function(data, t){
  sink(nullfile(), type="output")
  res = ComparisonSurv::Fixpoint.test(data$time, data$status, data$group, t0 = t)
  sink()
  return(res)
}

# Calculation of results
amll_summ = summary(amll_fit, times=20)
surv_KM = amll_summ$surv
se_KM = amll_summ$std.err
results = list(exp_surv = apply_all_tests(exp_surv, t = 1),
               amll = apply_all_tests(amll, t = 20))

# Calculation of reference results
ref_Comp = list( # ComparisonSurv reference
  exp_surv = apply_ComparisonSurv_tests(data=exp_surv, t=1),
  amll = apply_ComparisonSurv_tests(data=amll, t=20)
)

ref_bpcp = list( # bpcp reference
  exp_surv = apply_bpcp_tests(data=exp_surv, t=1),
  amll = apply_bpcp_tests(data=amll, t=20)
)

# The tests
testthat::test_that(
  "correct p-values and statistics for two numeric survival rates - data frame",
  {
    expect_equal(ref_Comp$exp_surv$test$pvalue[3:5],
                 results$exp_surv$p.values[3:5],
                 tolerance = 1e-5)
    expect_equal(ref_Comp$amll$test$statistic[3:5],
                 results$amll$statistics[3:5],
                 tolerance = 1e-5)

    expect_equal(ref_bpcp$exp_surv$p.values, results$exp_surv$p.values[1:3])
    expect_equal(ref_bpcp$amll$p.values, results$amll$p.values[1:3])
  }
)

testthat::test_that(
  "correct p-values and statistics for two numeric survival rates - survfit object",
  {
    results_fit = apply_all_tests(data=amll_fit, t = 20)
    expect_equal(ref_Comp$amll$test$statistic[3:5],
                 results_fit$statistics[3:5],
                 tolerance = 1e-5)

    expect_equal(ref_bpcp$exp_surv$p.values, results$exp_surv$p.values[1:3])
    expect_equal(ref_bpcp$amll$p.values, results$amll$p.values[1:3])
  }
)



testthat::test_that(
  "correct p-values and statistics for two numeric survival rates",
  {
    results_2num=apply_all_tests(data=NULL,surv_KM=surv_KM,se_KM=se_KM)
    expect_equal(ref_Comp$amll$test$statistic[3:5],
                 results_2num$statistics[3:5],
                 tolerance = 1e-5)
    expect_equal(ref_Comp$amll$test$pvalue[3:5],
                 results_2num$p.values[3:5],
                 tolerance = 1e-5)
    expect_equal(ref_bpcp$amll$p.values, results_2num$p.values[1:3])
  })

testthat::test_that(
  "correct p-values and statistics for one-group data frame and one-group numeric survival rate",
  {
    amll %>% filter(group==1) -> amll1
    results_1num=apply_all_tests(data=amll1,surv_KM=surv_KM[2],se_KM=se_KM[2],t=20)
    expect_equal(ref_Comp$amll$test$statistic[3:5],
                 results_1num$statistics[3:5],
                 tolerance = 1e-5)
    expect_equal(ref_Comp$amll$test$pvalue[3:5],
                 results_1num$p.values[3:5],
                 tolerance = 1e-5)
    expect_equal(ref_bpcp$amll$p.values, results_1num$p.values[1:3])
  })

testthat::test_that(
  "correct output of data frame names",
  {
    expect_equal(results$exp_surv$data.names, rep("data", 5))
    direct_data.names = c(naive.test(exp_surv, t = 1)$data.name,
                          logtra.test(exp_surv, t = 1)$data.name,
                          clog.test(exp_surv, t = 1)$data.name,
                          asinsqrt.test(exp_surv, t = 1)$data.name,
                          logit.test(exp_surv, t = 1)$data.name)
    expect_equal(direct_data.names, rep("exp_surv", 5))
  })

testthat::test_that(
  "correct output for individually defined variable names for time, group and status",
  {
    data = rotterdam
    t = 3000
    naive = naive.test(data, t = t, time = dtime, group = chemo, status = death)
    logtra = logtra.test(data, t = t, time = dtime, group = chemo, status = death)
    clog = clog.test(data, t = t, time = dtime, group = chemo, status = death)
    asinsqrt = asinsqrt.test(data, t = t, time = dtime, group = chemo, status = death)
    logit = logit.test(data, t = t, time = dtime, group = chemo, status = death)
    p.values = c(naive$p.value,
                 logtra$p.value,
                 clog$p.value,
                 asinsqrt$p.value,
                 logit$p.value)
    statistics = c(
      naive$statistic,
      logtra$statistic,
      clog$statistic,
      asinsqrt$statistic,
      logit$statistic
    )
    names(p.values) = NULL
    names(statistics) = NULL
    results_rttrdm = list(p.values = p.values, statistics = statistics)
    data %>% rename(group = chemo, status = death, time = dtime) -> data_renamed
    # ComparisonSurv reference
    ref_Comp_rttrdm = apply_ComparisonSurv_tests(data=data_renamed, t=t)
    # bpcp reference
    ref_bpcp_rttrdm = apply_bpcp_tests(data=data_renamed, t=t)
    expect_equal(ref_Comp_rttrdm$test$statistic[3:5],
                 results_rttrdm$statistics[3:5],
                 tolerance = 1e-5)
    expect_equal(ref_Comp_rttrdm$test$pvalue[3:5],
                 results_rttrdm$p.values[3:5],
                 tolerance = 1e-5)
    expect_equal(ref_bpcp_rttrdm$p.values, results_rttrdm$p.values[1:3])
  })

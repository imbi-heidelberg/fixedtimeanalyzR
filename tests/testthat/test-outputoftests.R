# In this statistical test, we test the output of the statistical tests.
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

# Preparations of data franes
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
  naive = fixtdiff(
    data$time,
    data$status,
    data$group,
    trans = "identity",
    varpooled = FALSE,
    testtime = t
  )
  logtra = fixtdiff(
    data$time,
    data$status,
    data$group,
    trans = "log",
    varpooled = FALSE,
    testtime = t
  )
  clog = fixtdiff(
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

# Calculation of results
amll_summ = summary(amll_fit, times=20)
surv_KM = amll_summ$surv
se_KM = amll_summ$std.err
results = list(exp_surv = apply_all_tests(exp_surv, t = 1),
               amll = apply_all_tests(amll, t = 20))

# Calculation of reference results
sink(nullfile())
ref_Comp = list( # ComparisonSurv reference
  exp_surv = Fixpoint.test(exp_surv$time, exp_surv$status, exp_surv$group, t0 =
                             1),
  amll = Fixpoint.test(amll$time, amll$status, amll$group, t0 =
                         20)
)
sink()

ref_bpcp = list( # bpcp reference
  exp_surv = apply_bpcp_tests(data=exp_surv, t=1),
  amll = apply_bpcp_tests(data=amll, t=20)
)

# The tests
test_that(
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

test_that(
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



test_that(
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

test_that(
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

test_that(
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

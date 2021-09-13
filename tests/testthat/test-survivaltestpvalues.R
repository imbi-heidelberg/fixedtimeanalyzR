context("survival test p-values")

library(dplyr)
data(exp_surv)
# ComparisonSurv only accepts 0 and 1 as group variables
exp_surv %>% mutate(group = group-1) -> exp_surv
naive = naive.test(exp_surv, t=1)
logtra = logtra.test(exp_surv, t=1)
clog = clog.test(exp_surv, t=1)
asinsqrt = asinsqrt.test(exp_surv,t=1)
logit = logit.test(exp_surv,t=1)
p.values = c(naive$p.value, logtra$p.value, clog$p.value, asinsqrt$p.value, logit$p.value)
statistics = c(naive$statistic, logtra$statistic, clog$statistic, asinsqrt$statistic, logit$statistic)
# Remove names, as references use different names
names(p.values) = NULL
names(statistics) = NULL

test_that("the statistics and p-values for clog.test, asinsqrt.test and
          logit.test are equal to the corresponding tests from the ComparisonSurv package", {
  library(ComparisonSurv)
            # TODO! Suppress the output of Fixpoint.test()
  reference = Fixpoint.test(exp_surv$time, exp_surv$status, exp_surv$group, t0=1)
  reference_p.values = reference$test[["pvalue"]]
  reference_statistics = reference$test[["statistic"]]
  # ComparisonSurv has two issues: (1) apparently, the statistics for naive and logtra were switched, so we omit these
  # (2) the p-values are only given up to the fourth digit.
  expect_equal(reference_p.values[3:5], p.values[3:5], tolerance=1e-5)
  expect_equal(reference_statistics[3:5], statistics[3:5], tolerance=1e-5)
})

test_that("the p-values for naive.test, logtra.test and clog.test are equal to the corresponding tests from bpcp package", {
  library(bpcp)

  reference_naive = fixtdiff(exp_surv$time, exp_surv$status, exp_surv$group, trans="identity", varpooled=FALSE, testtime=1)
  reference_logtra = fixtdiff(exp_surv$time, exp_surv$status, exp_surv$group, trans="log", varpooled=FALSE, testtime=1)
  reference_clog = fixtdiff(exp_surv$time, exp_surv$status, exp_surv$group, trans="cloglog", varpooled=FALSE, testtime=1)
  reference_p.values = c(reference_naive$p2, reference_logtra$p2, reference_clog$p2)
  names(reference_p.values) = NULL
  expect_equal(reference_p.values, p.values[1:3])
})

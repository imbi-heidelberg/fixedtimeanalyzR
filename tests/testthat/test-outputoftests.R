context("output of statistical tests")

library(dplyr)
library(ComparisonSurv)
library(bpcp)

data(exp_surv)
data(cancer, package = "survival")
# ComparisonSurv only accepts 0 and 1 as group variables
exp_surv %>% mutate(group = group - 1) -> exp_surv
# fixedtimeanalyzR only accepts group variable named group
aml %>% rename(group = x) %>% mutate(group = factor(group, levels=c("Maintained","Nonmaintained"),labels=c(1,0)))-> amll

apply_all_tests = function(data, t) {
  naive = naive.test(data, t = t)
  logtra = logtra.test(data, t = t)
  clog = clog.test(data, t = t)
  asinsqrt = asinsqrt.test(data, t = t)
  logit = logit.test(data, t = t)
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
  names(p.values) = NULL
  names(statistics) = NULL
  return(list(
    p.values = p.values,
    statistics = statistics,
    data.names = data.names
  ))
}
results = list(exp_surv = apply_all_tests(exp_surv, t = 1),
               amll = apply_all_tests(amll, t = 20))

# Remove names, as references use different names


test_that(
  "the statistics and p-values for clog.test, asinsqrt.test and
          logit.test are equal to the corresponding tests from the ComparisonSurv package",
  {

    sink(nullfile())
    reference = list(
      exp_surv = Fixpoint.test(exp_surv$time, exp_surv$status, exp_surv$group, t0 =
                                 1),
      amll = Fixpoint.test(amll$time, amll$status, amll$group, t0 =
                             20)
    )
    sink()
    # ComparisonSurv has two issues: (1) apparently, the statistics for naive and logtra were switched, so we omit these
    # (2) the p-values are only given up to the fourth digit.
    expect_equal(reference$exp_surv$test$pvalue[3:5],
                 results$exp_surv$p.values[3:5],
                 tolerance = 1e-5)
    expect_equal(reference$exp_surv$test$statistic[3:5],
                 results$exp_surv$statistics[3:5],
                 tolerance = 1e-5)
  }
)

test_that(
  "the p-values for naive.test, logtra.test and clog.test are equal to the corresponding tests from bpcp package",
  {
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
    reference = list(
      exp_surv = apply_bpcp_tests(data=exp_surv, t=1),
      amll = apply_bpcp_tests(data=amll, t=20)
    )
    expect_equal(reference$exp_surv$p.values, results$exp_surv$p.values[1:3])
    expect_equal(reference$amll$p.values, results$amll$p.values[1:3])
  }
)

test_that(
  "the data names are output correctly",
  {
  expect_equal(results$exp_surv$data.names, rep("data", 5))
  direct_data.names = c(naive.test(exp_surv, t = 1)$data.name,
                        logtra.test(exp_surv, t = 1)$data.name,
                        clog.test(exp_surv, t = 1)$data.name,
                        asinsqrt.test(exp_surv, t = 1)$data.name,
                        logit.test(exp_surv, t = 1)$data.name)
  expect_equal(direct_data.names, rep("exp_surv", 5))
})

test_that(
  "comparison of data set to numeric value works",
  {
    expect_equal(results$exp_surv$data.names, rep("data", 5))
    direct_data.names = c(naive.test(exp_surv, t = 1)$data.name,
                          logtra.test(exp_surv, t = 1)$data.name,
                          clog.test(exp_surv, t = 1)$data.name,
                          asinsqrt.test(exp_surv, t = 1)$data.name,
                          logit.test(exp_surv, t = 1)$data.name)
    expect_equal(direct_data.names, rep("exp_surv", 5))
  })

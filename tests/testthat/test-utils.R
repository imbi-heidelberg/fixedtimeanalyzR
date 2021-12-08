context("helper functions in survival_utils.R")

testthat::test_that(
  "correct output for individually defined variable names for time, group and status",
  {
    library(survival)
    data(cancer)
    # Existing group variable
    expect_equal(toString(get_survfit(data=rotterdam, time=dtime, status=death, group=chemo)$call),
                 "survfit, survival::Surv(dtime, death) ~ chemo, data")
    # Non-existant variables
    expect_error(get_survfit(data=rotterdam, time=dtime, status=death, group=nonexist))
    expect_error(get_survfit(data=rotterdam, time=dtime, status=nonexist, group=chemo))
    expect_error(get_survfit(data=rotterdam, time=nonexist, status=death, group=chemo))
    # Default values for time and status
    expect_equal(toString(get_survfit(data=aml,group=x)$call),
                 "survfit, survival::Surv(time, status) ~ x, data")
    # Default values for group! (TODO!)
    expect_equal(toString(get_survfit(data=aml, group=1)$call),
                 "survfit, survival::Surv(time, status) ~ 1, data")
    expect_warning(expect_equal(toString(get_survfit(data=aml, group=NULL)$call),
                 "survfit, survival::Surv(time, status) ~ 1, data"))
    # Convert time = NULL to default value
    expect_warning(expect_equal(toString(get_survfit(data=aml, time=NULL, status=status, group=x)$call),
                 "survfit, survival::Surv(time, status) ~ x, data"))
    # Convert status = NULL to the default value
    expect_warning(expect_equal(toString(get_survfit(data=aml, time=time, status=NULL, group=x)$call),
                 "survfit, survival::Surv(time, status) ~ x, data"))
  })

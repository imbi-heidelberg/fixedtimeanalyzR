context("helper function in survival_utils.R")

testthat::test_that(
  "correct output for individually defined variable names for time, group and status",
  {
    library(survival)
    data(cancer)
    # Existing group variable
    expect_equal(get_survfit(data=rotterdam, time=dtime, status=death, group=chemo),
                 "survfit, survival::Surv(dtime, death) ~ chemo, data")
    # Non-existant variables
    expect_error(get_survfit(data=rotterdam, time=dtime, status=death, group=nonexist))
    expect_error(get_survfit(data=rotterdam, time=dtime, status=nonexist, group=chemo))
    expect_error(get_survfit(data=rotterdam, time=nonexist, status=death, group=chemo))
    # Default values for time and status
    expect_equal(get_survfit(data=aml,group=x),
                 "survfit, survival::Surv(time, status) ~ x, data")
    # Default values for group! (TODO!)
    expect_equal(get_survfit(data=aml, group=1),
                 "survfit, survival::Surv(time, status) ~ 1, data")
    expect_equal(get_survfit(data=aml, group=NULL),
                 "survfit, survival::Surv(time, status) ~ 1, data")
    # Convert time = NULL to default value
    expect_equal(get_survfit(data=aml, time=NULL, status=status, group=x),
                 "survfit, survival::Surv(time, status) ~ x, data")
    # Convert status = NULL to the default value
    expect_equal(get_survfit(data=aml, time=time, status=NULL, group=x),
                 "survfit, survival::Surv(time, status) ~ x, data")
  })

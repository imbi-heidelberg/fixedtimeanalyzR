context("helper function in survival_utils.R")

testthat::test_that(
  "correct output for individually defined variable names for time, group and status",
  {
    library(survival)
    data(cancer)
    # Existing group variable
    get_survfit(data=rotterdam, time=dtime, status=death, group=chemo)
    # Non-existant group variable
    get_survfit(data=rotterdam, time=dtime, status=death, group=nonexist)
    # Default values for time and status
    get_survfit(data=aml,group=x)
    # Default values for group! (TODO!)

    # time or status = NULL
    get_survfit(data=rotterdam, time=NULL, status=death, group=nonexist)
  })

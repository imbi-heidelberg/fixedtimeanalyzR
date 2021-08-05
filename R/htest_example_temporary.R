farrington.manning <- function(
  group1,
  group2,
  delta = 0, # note that the null is formulates in terms of -delta!
  alternative = "greater",
  alpha = 0.025
) {
  # remove na's
  if (any(is.na(group1)) | any(is.na(group2))) {
    warning("There are missing values in either of the groups, will be ignored.")
    group1 <- group1[!is.na(group1)]
    group2 <- group2[!is.na(group2)]
  }



  # check for logical input vectors
  if (!is.logical(group1) | !is.logical(group2)) {
    stop("Inputs group1 and group2 must be logical vectors!")
  }



  # check alternative
  if (!(alternative %in% c("two.sided", "less", "greater"))) {
    stop("alternative must be either two.sided, less or greater")
  }



  # check alpha range
  if (0 >= alpha | alpha >= 1) {
    stop("alpha must be in (0, 1)")
  }



  # check delta
  if (delta >= 1) {
    stop("A rate difference of < 1 is not possible.")
  }
  if (delta <= -1) {
    stop("A rate difference of > 1 is not possible.")
  }



  # construct results object
  res <- list()
  class(res) <- "htest"
  res$null.value  <- c(delta)
  names(res$null.value) <- c("rate difference (group 1 - group 2)")
  res$alternative <- alternative
  str <- "Farrington-Manning test for rate differences"
  if (alternative == "greater" & delta < 0) {
    str <- "Non-inferiority test for rates according to Farrington-Manning"
  }
  if (alternative == "greater" & delta >= 0) {
    str <- "Superiorty test for rates according to Farrington-Manning"
  }
  res$method      <- "Farrington-Manning test for non-inferiority of rates"
  res$data.name   <- sprintf("group 1: %s, group 2: %s", deparse(substitute(group1)), deparse(substitute(group2)))
  res$parameters  <- c(delta)
  names(res$parameters) <- c("noninferiority margin")



  # number of samples per group
  n1 <- length(group1)
  n2 <- length(group2)
  res$sample.size <- n1 + n2



  # compute maximum likelihood estimates
  p1_ML <- mean(group1)
  p2_ML <- mean(group2)



  # rate difference
  diff_ML <- p1_ML - p2_ML
  res$estimate <- c(diff_ML)
  names(res$estimate) <- c("rate difference (group 1 - group 2)")



  # standard deviation of the rate difference under the null hypothesis (risk difference = -delta)
  get_sd_diff_ML_null <- function(delta) {
    theta           <- n2/n1
    d               <- -p1_ML*delta*(1 + delta)
    c               <- delta^2 + delta*(2*p1_ML + theta + 1) + p1_ML + theta*p2_ML
    b               <- -(1 + theta + p1_ML + theta*p2_ML + delta*(theta + 2))
    a               <- 1 + theta
    v               <- b^3/(27*a^3) - b*c/(6*a^2) + d/(2*a)
    u               <- sign(v)*sqrt(b^2/(9*a^2) - c/(3*a))
    w               <- (pi + acos(v/u^3))/3
    p1_ML_null      <- 2*u*cos(w) - b/(3*a)
    p2_ML_null      <- p1_ML_null - delta
    sd_diff_ML_null <- sqrt(p1_ML_null*(1 - p1_ML_null)/n1 + p2_ML_null*(1 - p2_ML_null)/n2)
    return(sd_diff_ML_null)
  }
  sd_diff_ML_null <- get_sd_diff_ML_null(delta)



  # test statistic
  get_z <- function(delta) {
    z <- (diff_ML - delta)/get_sd_diff_ML_null(delta)
    return(z)
  }
  z <- get_z(delta)
  res$statistic <- z
  names(res$statistic) <- "Z-statistic"



  # p-value, probability of Z > < == z
  p_value_greater <- 1 - pnorm(z)
  p_value_less <- pnorm(z)
  p_value_two.sided <- 2*min(p_value_less, p_value_greater)
  if (alternative == "greater") {
    res$p.value <- p_value_greater
  }
  if (alternative == "less") {
    res$p.value <- p_value_less
  }
  if (alternative == "two.sided") {
    res$p.value <- p_value_two.sided
  }



  # confidence interval by inversion of two-sided test
  p_value_two.sided <- function(delta) {
    z <- get_z(delta)
    p_value_greater <- 1 - pnorm(z)
    p_value_less <- pnorm(z)
    2*min(p_value_less, p_value_greater)
  }



  alpha_mod <- ifelse(alternative == "two.sided", alpha, 2*alpha)
  ci_lo <- uniroot(
    function(delta) p_value_two.sided(delta) - alpha_mod, interval = c(-1+1e-6, res$estimate), tol = 1e-12
  )$root
  ci_hi <- uniroot(
    function(delta) p_value_two.sided(delta) - alpha_mod, interval = c(res$estimate, 1 - 1e-6), tol = 1e-12
  )$root



  # confidence interval
  res$conf.int <- c(ci_lo, ci_hi)
  attr(res$conf.int, "conf.level") <- 1 - alpha_mod



  return(res)
}

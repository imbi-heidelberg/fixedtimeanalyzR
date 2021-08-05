naive.test <- function(data, surv_KM, se_KM, t) {
  # NA warning.
  if (any(is.na(data))) {
    warning("There are missing values in either of the groups.")
    #data <- data[!is.na(data)]
  }

  # TODO: Im Test sollte es drei verschiedene Möglichkeiten geben. Entweder ein
  # Data.frame mit zwei Gruppen, oder ein Data.frame mit nur einer Gruppe +
  # surv_KM und se_KM mit jeweils einer Komponente, oder kein Data.frame und
  # surv_KM und se_KM mit jeweils zwei Komponenten.

  # construct results object
  res <- list()
  class(res) <- "htest"
  res$null.value  <- c(0)
  names(res$null.value) <- c("absolute survival rate difference |group 1 - group 2|")
  res$alternative <- "greater" # A chi-squared test is always a one-sided test
  res$method <- "Naive test for equality of survival rates at a fixed point in time.

  This test is using the naive test statistic X1^2 from Klein et al.'s paper
  'Analyzing survival curves at a fixed point in time', published in Stat. Med.,
  26 (2007)."
  res$data.name   <- sprintf(deparse(substitute(data)))
  res$parameters  <- c(t)
  names(res$parameters) <- c("time")

  # Survival rates and their difference
  surv_KM <- get_surv_KM(data, t)
  diff_surv_KM <- abs(surv_KM[1] - surv_KM[2])
  res$estimate <- c(surv_KM[1], surv_KM[2], diff_surv_KM)
  names(res$estimate) <- c("Survival rate of group 1", "Survival rate of group 2", "Absolute difference")

  # Test statistic
  t <- naive.t(data, t)
  res$statistic <- t
  names(res$statistic) <- "naive statistic X1^2"

  # p-value, probability of naive.t < t
  res$p.value <- 1 - pchisq(t, df = 1)
  return(res)
}

# TODO: Was ist, wenn der NutzerIn keine Daten, sondern nur Überlebensraten eingeben möchte.

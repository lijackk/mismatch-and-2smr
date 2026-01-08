#Simex as applied to regression
simex_regression <- function(predictor, predictor.se, outcome, outcome.se, lambda_vec = seq(0, 5, by = 0.25), iter = 1000) {
  stopifnot(length(predictor) == length(predictor.se) & length(predictor) == length(outcome) & length(predictor) == length(outcome.se))
  simex_estimates <- vector()
  N <- length(predictor)
  for (i in 1:iter) {
    slope_lambda <- vector()
    for (j in 1:length(lambda_vec)) {
      bootstrap_resample <- sample(N, N, replace = TRUE)
      lambda <- lambda_vec[j]
      new_predictor <- rnorm(predictor[bootstrap_resample],
                         mean = predictor[bootstrap_resample],
                         sd = sqrt(1 + lambda)*predictor.se[bootstrap_resample])
      new_outcome <- rnorm(outcome[bootstrap_resample],
                         mean = outcome[bootstrap_resample],
                         sd = sqrt(1 + lambda)*outcome.se[bootstrap_resample])
      slope_lambda[j] <- lm(new_outcome ~ new_predictor - 1, weights = 1/outcome.se[bootstrap_resample]^2)$coefficients
    }
    simex_model <- lm(slope_lambda ~ lambda_vec)$coefficients
    simex_estimates[i] <- simex_model[1] + (-1)*simex_model[2]
    if (i %% 100 == 0) print(i)
  }
  return(simex_estimates)
}

#Simex as applied to regression with an intercept term
simex_regression_intercept <- function(predictor, predictor.se, outcome, outcome.se, lambda_vec = seq(0, 5, by = 0.25), iter = 1000) {
  stopifnot(length(predictor) == length(predictor.se) & length(predictor) == length(outcome) & length(predictor) == length(outcome.se))
  simex_slope_estimates <- vector()
  simex_intercept_estimates <- vector()
  N <- length(predictor)
  for (i in 1:iter) {
    slope_lambda <- vector()
    intercept_lambda <- vector()
    for (j in 1:length(lambda_vec)) {
      bootstrap_resample <- sample(N, N, replace = TRUE)
      lambda <- lambda_vec[j]
      new_predictor <- rnorm(predictor[bootstrap_resample],
                         mean = predictor[bootstrap_resample],
                         sd = sqrt(1 + lambda)*predictor.se[bootstrap_resample])
      new_outcome <- rnorm(outcome[bootstrap_resample],
                         mean = outcome[bootstrap_resample],
                         sd = sqrt(1 + lambda)*outcome.se[bootstrap_resample])
      slope_lambda[j] <- lm(new_outcome ~ new_predictor, weights = 1/outcome.se[bootstrap_resample]^2)$coefficients[2]
      intercept_lambda[j] <- lm(new_outcome ~ new_predictor, weights = 1/outcome.se[bootstrap_resample]^2)$coefficients[1]
    }
    simex_slope_model <- lm(slope_lambda ~ lambda_vec)$coefficients
    simex_slope_estimates[i] <- simex_slope_model[1] + (-1)*simex_slope_model[2]
  
    simex_intercept_model <- lm(intercept_lambda ~ lambda_vec)$coefficients
    simex_intercept_estimates[i] <- simex_intercept_model[1] + (-1)*simex_intercept_model[2]
    if (i %% 100 == 0) print(i)
  }
  return(list(slope = simex_slope_estimates, intercept = simex_intercept_estimates))
}

#Simex regression, one iteration, no bootstrapping
  #Returns the point estimate, and lambda vector
simex_regression.nb <- function(predictor, predictor.se, outcome, outcome.se, lambda_vec = seq(0, 5, by = 0.25)) {
  stopifnot(length(predictor) == length(predictor.se) & length(predictor) == length(outcome) & length(predictor) == length(outcome.se))
  N <- length(predictor)
  slope_lambda <- vector()
  for (j in 1:length(lambda_vec)) {
      lambda <- lambda_vec[j]
      new_predictor <- rnorm(predictor,
                         mean = predictor,
                         sd = sqrt(1 + lambda)*predictor.se)
      new_outcome <- rnorm(outcome,
                         mean = outcome,
                         sd = sqrt(1 + lambda)*outcome.se)
      slope_lambda[j] <- lm(new_outcome ~ new_predictor - 1, weights = 1/outcome.se^2)$coefficients
    }
  simex_model <- lm(slope_lambda ~ lambda_vec)$coefficients
  simex_estimate <- simex_model[1] + (-1)*simex_model[2]
  return(list(simex_estimate = simex_estimate, slope_lambda = slope_lambda, lambda_vec = lambda_vec, simex_model = simex_model))
}


#Simex for a median with measurement error
simex_median <- function(num, num.se, denom, denom.se, lambda_vec = seq(0, 5, by = 0.25), iter = 1000) {
  stopifnot(length(num) == length(num.se) & length(num) == length(denom) & length(num) == length(denom.se))
  simex_estimates <- vector()
  N <- length(num)
  for (i in 1:iter) {
    median_lambda <- vector()
    for (j in 1:length(lambda_vec)) {
      bootstrap_resample <- sample(N, N, replace = TRUE)
      lambda <- lambda_vec[j]
      new_num <- rnorm(num[bootstrap_resample],
                         mean = num[bootstrap_resample],
                         sd = sqrt(1 + lambda)*num.se[bootstrap_resample])
      new_denom <- rnorm(denom[bootstrap_resample],
                         mean = denom[bootstrap_resample],
                         sd = sqrt(1 + lambda)*denom.se[bootstrap_resample])
      median_lambda[j] <- median(new_num/new_denom)
    }
    simex_model <- lm(median_lambda ~ lambda_vec)$coefficients
    simex_estimates[i] <- simex_model[1] + (-1)*simex_model[2]
    if (i %% 100 == 0) print(i)
  }
  return(simex_estimates)
}

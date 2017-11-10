context("corrected posterior")


test_that("'compute_ratio_priors_gaussian' works as expected", {

  # should be ok
  x= seq(-6, 6, length.out=10)
  inla_post = matrix(c(x, dnorm(x)), 10, 2)
  x = compute_ratio_priors_gaussian(hyperpar_new = c(mean = 1, prec = 1), 
                                    hyperpar_old = c(mean = 1, prec = 1.1), 
                                    inla_internal_marginal = inla_post)
  x_test = c(11.0490551, 4.7488404, 2.4381360, 1.4953270, 1.0955251, 0.9587743, 
             1.0023477, 1.2517807, 1.8674374, 3.3279114)
  expect_equal(x, x_test, tolerance = 1e-3)
  
})

test_that("'compute_ratio_priors_gamma' works as expected", {
  
# should be ok
x = seq(-6, 6, length.out=10)
inla_post = matrix(c(x, dnorm(x)), 10, 2)
x = compute_ratio_priors_gamma(hyperpar_new = c(shape = 1, rate = 1),  
                               hyperpar_old = c(shape = 1, rate = 1.01),  
                               inla_internal_marginal = inla_post)
x_test = c(0.9901236, 0.9901921, 0.9904523,  0.9914399,  0.9951954,  1.0095725,  
           1.0660287,  1.3104439,  2.8676028, 55.9432399)  
expect_equal(x, x_test, tolerance = 1e-3)

})

test_that("'compute_ratio_priors_custom' works as expected", {
  
  # user defined gamma log prior density
  log_f <- function(x, hyperpar){
    dgamma(exp(x), shape = exp(hyperpar[1]), rate = exp(hyperpar[2]), log = TRUE) + x
  }
  
  # should be ok
  x = seq(-6, 6, length.out=10)
  inla_post = matrix(c(x, dnorm(x)), 10, 2)
  x = compute_ratio_priors_custom(hyperpar_new = c(log(1), log(1)),  
                                 hyperpar_old = c(log(1), rate = log(1.01)),  
                                 inla_internal_marginal = inla_post,
                                  log_prior_density = log_f)
  
  x_test = c(0.9901236, 0.9901921, 0.9904523,  0.9914399,  0.9951954,  1.0095725,  
             1.0660287,  1.3104439,  2.8676028, 55.9432399)  
  expect_equal(x, x_test, tolerance = 1e-3)
  
})

test_that("'compute_corrected_posteriors' works as expected", {

  x = seq(-6, 6, length.out=10)
  inla_post = matrix(c(x, dnorm(x)), 10, 2)
  
  ## test 'prior_type' argument
  
  # should be error - "compute_corrected_posteriors: Currently 'prior_type' needs to be 'gaussian' or 'gamma'."
  expect_error(compute_corrected_posteriors(prior_type = NULL,
                               hyperpar_new = c(mean = 1, prec = 1.01),
                               hyperpar_old = c(mean = 1, prec = 1),
                               inla_marginal_posterior = inla_post, 
                               integration_limits = TRUE),
               "compute_corrected_posteriors: Currently 'prior_type' needs to be 'gaussian', 'gamma' or 'custom'.")
  
  # should be error - "compute_corrected_posteriors: Currently 'prior_type' needs to be 'gaussian' or 'gamma'."
  expect_error(compute_corrected_posteriors(prior_type = "CRAZY_INPUT",
                               hyperpar_new = c(mean = 1, prec = 1.01),
                               hyperpar_old = c(mean = 1, prec = 1),
                               inla_marginal_posterior = inla_post, 
                               integration_limits = TRUE),
               "compute_corrected_posteriors: Currently 'prior_type' needs to be 'gaussian', 'gamma' or 'custom'.")
  
  ## test gaussian case arguments
  
  # should be ok
  x = compute_corrected_posteriors(prior_type = "gaussian",
                                   hyperpar_new = c(mean = 1, prec = 1.01),
                                   hyperpar_old = c(mean = 1, prec = 1),
                                   inla_marginal_posterior = inla_post, 
                                   integration_limits = TRUE)
  expect_true(ifelse(is.matrix(x), ncol(x) == 2, FALSE))
  
  ## test gamma case arguments
  
  # should be ok
  x = compute_corrected_posteriors(prior_type = "gamma",
                                   hyperpar_new = c(shape = 1, rate = 1.01),
                                   hyperpar_old = c(shape = 1, rate = 1),
                                   inla_marginal_posterior = inla_post, 
                                   integration_limits = TRUE)
  expect_true(ifelse(is.matrix(x), ncol(x) == 2, FALSE))

  ## test custom case arguments

  # user defined gamma log prior density
  log_f <- function(x, hyperpar){
    dgamma(exp(x), shape = exp(hyperpar[1]), rate = exp(hyperpar[2]), log = TRUE) + x
  }
  
  # should be ok
  x = compute_corrected_posteriors(prior_type = "custom",
                                   hyperpar_new = c(log(1), log(1.01)),
                                   hyperpar_old = c(log(1), log(1)),
                                   inla_marginal_posterior = inla_post, 
                                   integration_limits = TRUE,
                                   log_prior_density = log_f)
  expect_true(ifelse(is.matrix(x), ncol(x) == 2, FALSE))
  
  
  ## test integration_limits = FALSE
  
  # should be ok
  x = compute_corrected_posteriors(prior_type = "gaussian",
                                   hyperpar_new = c(mean = 1, prec = 1.01),
                                   hyperpar_old = c(mean = 1, prec = 1),
                                   inla_marginal_posterior = inla_post, 
                                   integration_limits = FALSE)
  expect_true(ifelse(is.matrix(x), ncol(x) == 2, FALSE))

})

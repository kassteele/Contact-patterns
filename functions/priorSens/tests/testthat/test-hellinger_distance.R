context("Helliner distance")

# Test cases for Helliner distance
#----------------------------------

test_that("'compute_hellinger' works as expected", {

  # (Gaussian, Gaussian) case

  ## correct use goes ok.
  hyperpar1 = cbind(mean = 0, prec = 1)
  hyperpar2 = cbind(mean = 1, prec = 1)

  expect_identical(compute_hellinger(log_dist1 = "gaussian", hyperpar1 = hyperpar1, 
                                     log_dist2 = "gaussian", hyperpar2 = hyperpar2),
                   compute_hellinger_gaussians(hyperpar1 = hyperpar1, 
                                               hyperpar2 = hyperpar2))

  ## error in the hyperparameter argument stops.
  hyperpar1 = c(mean = 0, 1)
  hyperpar2 = c(mean = 1, prec = 1)
  
#   expect_error(compute_hellinger(log_dist1 = "gaussian", hyperpar1 = hyperpar1, 
#                                  log_dist2 = "gaussian", hyperpar2 = hyperpar2),
#                "hyperpar should be a named numeric array containing 'mean' and 'prec' parameters for the Gaussian distribution. 'prec' should be > 0.")

  # (Gamma, Gamma) case

  ## correct use goes ok.
  hyperpar1 = cbind(shape = 1, rate = 0.05)
  hyperpar2 = cbind(shape = 1, rate = 0.2)
  
  expect_identical(compute_hellinger(log_dist1 = "gamma", hyperpar1 = hyperpar1, 
                                     log_dist2 = "gamma", hyperpar2 = hyperpar2),
                   compute_hellinger_gammas(hyperpar1 = hyperpar1, 
                                           hyperpar2 = hyperpar2))
  
  ## error in the hyperparameter argument stops.
#   hyperpar1 = c(shape = 1, 0.05)
#   hyperpar2 = c(shape = 1, rate = 0.2)
#   
#   expect_error(compute_hellinger(log_dist1 = "gamma", hyperpar1 = hyperpar1, 
#                                  log_dist2 = "gamma", hyperpar2 = hyperpar2),
#                "hyperpar should be a named numeric array containing 'shape' and 'rate' parameters for the Gamma distribution. Both shape and rate should be > 0.")

  # (function, function) case
  
  log_exponential <- function(x, hyperpar){
    dexp(x, rate = hyperpar, log = TRUE)
  }

  expect_equal(compute_hellinger(log_dist1 = log_exponential, hyperpar1 = 1, 
                                 log_dist2 = log_exponential, hyperpar2 = 2, domain = c(0, Inf)),
               0.2391463, tolerance = 1e-3)
#     
#   # should be error - "'domain' is required for numerically computed HD."
#   expect_error(compute_hellinger(log_dist1 = log_exponential, hyperpar1 = 1, 
#                                  log_dist2 = log_exponential, hyperpar2 = 2, domain = NULL),
#                "'domain' is required for numerically computed HD.")
  
})

test_that("'compute_hellinger_distance_inla' works as expected", {

  x= seq(-6, 6, length.out=100)
  # Should be approximately zero
  inla_post = matrix(c(x, dnorm(x)), 100, 2)
  
  expect_equal(compute_hellinger_distance_inla(inla_posterior_1 = inla_post, 
                                  inla_posterior_2 = inla_post, 
                                  integration_limits = TRUE, 
                                  renormalize = TRUE), 
               0,
               tolerance = 1e-3)
  
  # different domain, should issue an error
  inla_post1 = inla_post
  inla_post1[1,1] = inla_post1[1,1] - 0.1
  
  expect_error(compute_hellinger_distance_inla(inla_posterior_1 = inla_post, 
                                               inla_posterior_2 = inla_post1, 
                                               integration_limits = TRUE, 
                                               renormalize = TRUE), 
               "compute_hellinger_distance_inla: Domains of the inla marginals should be the same.")
  
  
})
  

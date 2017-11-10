context("sensitivity")


test_that("'compute_ratio_priors_gaussian' works as expected", {

  x = seq(-6, 6, length.out=10)
  inla_post = matrix(c(x, dnorm(x)), 10, 2)
  
  grid = list(polar = as.matrix(data.frame(theta = c(-3.141593, 0.000000, 3.141593), r = c(1.00002, 1.00002, 1.00002))),
              cartesian = as.matrix(data.frame(mean = c(0.8584916, 1.1415084, 0.8584916), prec = c(1, 1, 1))),
              hyperpar = c(mean = 1, prec = 1),
              grid_epsilon = 0.05,
              prior_type = "gaussian")
  class(grid) <- "grid"
  
  # should be ok
  sens_obj = compute_sensitivity(grid = grid, 
                                 inla_marginal_posterior = inla_post,
                                 integration_limits = TRUE)
  x = sens_obj$sensitivity
  expect_true(ifelse(is.matrix(x), ((all(colnames(x) %in% c("theta", "r", "mean", "prec", "sensitivity"))) & 
                                      all(dim(x) == c(3, 5))), FALSE))
  expect_identical(grid, sens_obj$polar_grid)
  
  # should be ok
  sens_obj = compute_sensitivity(grid = grid, 
                                 inla_marginal_posterior = inla_post,
                                 integration_limits = FALSE)
  x = sens_obj$sensitivity
  expect_true(ifelse(is.matrix(x), ((all(colnames(x) %in% c("theta", "r", "mean", "prec", "sensitivity"))) & 
                                      all(dim(x) == c(3, 5))), FALSE))
  
  # should be error - "'grid' argument needs to be an object with class 'grid'.'"
  expect_error(compute_sensitivity(grid = "CRAZY_INPUT", 
                      inla_marginal_posterior = inla_post,
                      integration_limits = TRUE),
               "'grid' argument needs to be an object with class 'grid'.'")
  
  # should be error = "'inla_internal_marginal' is supposed to be a matrix."
  expect_error(compute_sensitivity(grid = grid, 
                      inla_marginal_posterior = c(1,2,3),
                      integration_limits = TRUE),
               "'inla_marginal_posterior' is supposed to be a matrix.")
  
})
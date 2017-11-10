context("Grid")

# Test cases for Helliner distance
#----------------------------------

test_that("'coord_correc_func' works as expected", {

  cc = coord_correc_func(theta = pi/2,
                    coordinate_correction = NULL)
  cc_test = c(1, 1)
  expect_identical(cc, cc_test)

  cc = coord_correc_func(theta = -3*pi/4,
                         coordinate_correction = list(xneg = 1, xpos = 2, yneg = 3, ypos = 4))
  cc_test = c(1, 3)
  expect_identical(cc, cc_test)

  cc = coord_correc_func(theta = -pi/4,
                         coordinate_correction = list(xneg = 1, xpos = 2, yneg = 3, ypos = 4))
  cc_test = c(2, 3)
  expect_identical(cc, cc_test)

  cc = coord_correc_func(theta = pi/4,
                         coordinate_correction = list(xneg = 1, xpos = 2, yneg = 3, ypos = 4))
  cc_test = c(2, 4)
  expect_identical(cc, cc_test)

  cc = coord_correc_func(theta = 3*pi/4,
                         coordinate_correction = list(xneg = 1, xpos = 2, yneg = 3, ypos = 4))
  cc_test = c(1, 4)
  expect_identical(cc, cc_test)
  
})

test_that("'objective_to_compute_roots_r_given_theta_gaussian' works as expected", {

  # Should return 0.8315286
  expect_equal(objective_to_compute_roots_r_given_theta_gaussian(x = 2, 
                                                    hyperpar = c(mean = 2, prec = 3), 
                                                    fixed_theta = pi/2,
                                                    grid_epsilon = 0.05,
                                                    coordinate_correction = NULL),
               0.8315286, tolerance = 1e-3)

  # Should return 0.9322646
  expect_equal(objective_to_compute_roots_r_given_theta_gaussian(x = 2, 
                                                    hyperpar = c(mean = 2, prec = 3), 
                                                    fixed_theta = pi/2,
                                                    grid_epsilon = 0.05,
                                                    coordinate_correction = list(xneg = 1, xpos = 2,
                                                                                 yneg = 1, ypos = 2)),
               0.9322646, tolerance = 1e-3)

})

test_that("'objective_to_compute_roots_r_given_theta_gamma' works as expected", {
  
  # Should return 0.9487648
  expect_equal(objective_to_compute_roots_r_given_theta_gamma(x = 2, 
                                                 hyperpar = c(shape = 2, rate = 3), 
                                                 fixed_theta = pi/2,
                                                 grid_epsilon = 0.05,
                                                 coordinate_correction = NULL),
               0.9487648, tolerance = 1e-3)
  
  # Should return 0.9499992
  expect_equal(objective_to_compute_roots_r_given_theta_gamma(x = 2, 
                                                 hyperpar = c(shape = 2, rate = 3), 
                                                 fixed_theta = pi/2,
                                                 grid_epsilon = 0.05,
                                                 coordinate_correction = list(xneg = 1, xpos = 2,
                                                                              yneg = 1, ypos = 2)),
               0.9499992, tolerance = 1e-3)
  
})

test_that("'compute_roots_r_given_thetas' works as expected", {

  # gaussian case - should return c(1.5707963, 0.2002922)
  expect_equal(compute_roots_r_given_thetas(new_theta = pi/2,
                                            objective_function = objective_to_compute_roots_r_given_theta_gaussian,
                                            hyperpar = c(mean = 1, prec = 2),
                                            grid_epsilon = 0.05,
                                            coordinate_correction = NULL,
                                            method = "nlminb"),
               c(1.5707963, 0.2002922), tolerance = 1e-3)
  
  # gamma case - should return c(1.5707963, 0.1415689)
  expect_equal(compute_roots_r_given_thetas(new_theta = pi/2,
                                            objective_function = objective_to_compute_roots_r_given_theta_gamma,
                                            hyperpar = c(shape = 1, rate = 2),
                                            grid_epsilon = 0.05,
                                            coordinate_correction = NULL,
                                            method = "nlminb"),
               c(1.5707963, 0.1415689), tolerance = 1e-3)
  
  # custom case - should return c(1.5707963 0.3041684)
  log_gamma <- function(x, hyperpar){
    dgamma(x, shape = exp(hyperpar[1]), rate = exp(hyperpar[2]), log = TRUE)
  }
  
  expect_equal(compute_roots_r_given_thetas(new_theta = pi/2,
                                            objective_function = objective_to_compute_roots_r_given_theta_custom,
                                            hyperpar = c(log(1), log(2)),
                                            grid_epsilon = 0.05,
                                            coordinate_correction = NULL,
                                            log_prior_density = log_gamma,
                                            method = "nlminb"),
               c(1.5707963, 0.1415689), tolerance = 1e-3)
  
  # should return an error, new_theta needs to be between [-pi, pi]
  expect_error(compute_roots_r_given_thetas(new_theta = -2*pi,
                                            objective_function = objective_to_compute_roots_r_given_theta_custom,
                                            hyperpar = c(log(1), log(2)),
                                            grid_epsilon = 0.05,
                                            coordinate_correction = NULL,
                                            log_prior_density = log_gamma,
                                            method = "nlminb"))
  
  # should return an error, new_theta needs to be between [-pi, pi]
  expect_error(compute_roots_r_given_thetas(new_theta = 2*pi,
                                            objective_function = objective_to_compute_roots_r_given_theta_custom,
                                            hyperpar = c(log(1), log(2)),
                                            grid_epsilon = 0.05,
                                            coordinate_correction = NULL,
                                            log_prior_density = log_gamma,
                                            method = "nlminb"))
  
  # should return an error, grid_epsilon needs to be > 0. "'grid_epsilon' needs to be greater than zero."
  expect_error(compute_roots_r_given_thetas(new_theta = pi/2,
                                            objective_function = objective_to_compute_roots_r_given_theta_custom,
                                            hyperpar = c(log(1), log(2)),
                                            grid_epsilon = -0.05,
                                            coordinate_correction = NULL,
                                            log_prior_density = log_gamma,
                                            method = "nlminb"),
               "'grid_epsilon' needs to be greater than zero.")
  
  # should return an error - "'coordinate_correction' needs to be a list with the following ellements 'xneg', 'xpos', 'yneg', 'ypos'"
  expect_error(compute_roots_r_given_thetas(new_theta = pi/2,
                                            objective_function = objective_to_compute_roots_r_given_theta_custom,
                                            hyperpar = c(log(1), log(2)),
                                            grid_epsilon = 0.05,
                                            coordinate_correction = list(xpos = 1, xneg = 2,
                                                                         ypos = 3, kneg = 4),
                                            log_prior_density = log_gamma,
                                            method = "nlminb"),
               "'coordinate_correction' needs to be a list with the following ellements 'xneg', 'xpos', 'yneg', 'ypos'")
  
})

test_that("'compute_grid_polar' works as expected", {
  
  # gaussian case
  grid_test = list(polar = as.matrix(data.frame(theta = c(-3.141593, -3.141593, -1.570796, -1.047198, 0.000000, 1.047198, 1.570796, 3.141593),
                                                r = c(1.0000000, 1.000000, 1.000000, 1.0109638, 1.000000, 0.9893087, 1.000000, 1.0000000))),
                   cartesian = as.matrix(data.frame(mean = c(1.918299, 1.918299, 2.000000, 2.041298, 2.081701, 2.040414, 2.000000, 1.918299),
                                                    prec = c(3.000000, 3.000000, 2.455475, 2.517465, 3.000000, 3.561622, 3.665279, 3.000000))),
                   hyperpar = c(2, 3),
                   grid_epsilon = 0.05,
                   prior_type = "gaussian",
                   log_prior_density = NULL)
  class(grid_test) <- "grid"  

  grid = compute_grid_polar(number_axis_points = 4,
                            log_prior_density = "gaussian",
                            hyperpar = c(mean = 2, prec = 3),
                            grid_epsilon = 0.05)
  expect_equal(grid, grid_test, tolerance = 1e-3, check.attributes = TRUE, check.names = TRUE)
  
  # gamma case
  grid_test = list(polar = as.matrix(data.frame(theta = c(-3.141593, -3.141593, -1.570796, -1.047198, 0.000000, 1.047198, 1.570796, 3.141593),
                                                r = c(0.9999996, 1.0000000, 1.0000000, 0.7688738, 1.0000000, 1.7220470, 1.0000000, 0.9999996))),
                   cartesian = as.matrix(data.frame(shape = c(0.894051, 0.894051, 1.000000, 1.042679, 1.114844, 1.098127, 1.000000, 0.894051),
                                                    rate = c(2.000000, 2.000000, 1.735991, 1.820082, 2.000000, 2.470138, 2.304160, 2.000000))),
                   hyperpar = c(shape = 1, rate = 2),
                   grid_epsilon = 0.05,
                   prior_type = "gamma",
                   log_prior_density = NULL)
  class(grid_test) <- "grid"
  
  grid = compute_grid_polar(number_axis_points = 4,
                            log_prior_density = "gamma",
                            hyperpar = c(shape = 1, rate = 2),
                            grid_epsilon = 0.05)
  expect_equal(grid, grid_test, tolerance = 1e-3, check.attributes = TRUE, check.names = TRUE)

  # custom case
  log_loggamma <- function(x, hyperpar){
    dgamma(exp(x), shape = exp(hyperpar[1]), rate = exp(hyperpar[2]), log = TRUE) + x
  }
  
  grid_test = list(polar = as.matrix(data.frame(theta = c(-3.141593, -3.141593, -1.570796, -1.047198, 0.000000, 1.047198, 1.570796, 3.141593),
                                                r = c(0.9999996, 1.0000000, 1.0000000, 0.7688738, 1.0000000, 1.7220470, 1.0000000, 0.9999996))),
                   cartesian = as.matrix(data.frame(new_hyper1 = c(-1.119924e-01, -1.119925e-01, 6.857560e-18, 4.179380e-02, 1.087143e-01, 
                                                                   9.360558e-02, 6.656832e-18, -1.119924e-01),
                                                    new_hyper2 = c(0.6931472, 0.6931472, 0.5515783, 0.5988815, 0.6931472, 0.9042741, 0.8347161, 0.6931472))),
                   hyperpar = c(log(1), log(2)),
                   grid_epsilon = 0.05,
                   prior_type = "custom",
                   log_prior_density = log_loggamma)
  class(grid_test) <- "grid"
  
  grid = compute_grid_polar(number_axis_points = 4,
                            log_prior_density = log_loggamma,
                            hyperpar = c(log(1), log(2)),
                            grid_epsilon = 0.05)
  expect_equal(grid, grid_test, tolerance = 1e-3, check.attributes = TRUE, check.names = TRUE)  

  # expect warning - warning issued because the hyperpar is misspecified. It is hard to throw an
  # error here because a error is possible and NA is returned when that happens, so we want
  # the algorithm to keep going even if some errors occur.
#   expect_warning(compute_grid_polar(number_axis_points = 4,
#                      objective_function = "gaussian",
#                      hyperpar = c(mean = 2, 3),
#                      grid_epsilon = 0.05,
#                      search_interval = c(-100, 40)))

  # should return an error, "'number_axis_points' needs to be > 0."
  expect_error(compute_grid_polar(number_axis_points = 0,
                                  log_prior_density = "gaussian",
                     hyperpar = c(mean = 2, prec = 3),
                     grid_epsilon = 0.05), "'number_axis_points' needs to be > 0.")
  
  # should return an error, "Invalid 'objective_function'."
  expect_error(compute_grid_polar(number_axis_points = 4,
                                  log_prior_density = "CRAZY_INPUT",
                     hyperpar = c(mean = 2, prec = 3),
                     grid_epsilon = 0.05), "Invalid 'log_prior_density'.")
  
  # should return an error, "'grid_epsilon' needs to be greater than zero."
  expect_error(compute_grid_polar(number_axis_points = 4,
                                  log_prior_density = "gaussian",
                     hyperpar = c(mean = 2, prec = 3),
                     grid_epsilon = -0.05), "'grid_epsilon' needs to be greater than zero.")
  
})


context("Auxiliary functions")

test_that("'check_hyperpar_gaussian' works as expected", {

  hyperpar = c(mean = 1, prec = 1)
  expect_identical(check_hyperpar_gaussian(hyperpar), "ok")
  
  hyperpar = c(mean = 1, prec = 1, 3) 
  expect_identical(check_hyperpar_gaussian(hyperpar), "hyperpar should be a named numeric array containing 'mean' and 'prec' parameters for the Gaussian distribution. 'prec' should be > 0.")

  hyperpar = c(mean = 1, 1) 
  expect_identical(check_hyperpar_gaussian(hyperpar), "hyperpar should be a named numeric array containing 'mean' and 'prec' parameters for the Gaussian distribution. 'prec' should be > 0.")

  hyperpar = c(1, 1) 
  expect_identical(check_hyperpar_gaussian(hyperpar), "hyperpar should be a named numeric array containing 'mean' and 'prec' parameters for the Gaussian distribution. 'prec' should be > 0.")

  hyperpar = c(mean = 1, prec = -1)
  expect_identical(check_hyperpar_gaussian(hyperpar), "hyperpar should be a named numeric array containing 'mean' and 'prec' parameters for the Gaussian distribution. 'prec' should be > 0.")
  
})

test_that("'check_hyperpar_gamma' works as expected", {

  hyperpar = c(shape = 1, rate = 1)
  expect_identical(check_hyperpar_gamma(hyperpar), "ok")
  
  hyperpar = c(shape = 1, rate = 1, 3) 
  expect_identical(check_hyperpar_gamma(hyperpar), "hyperpar should be a named numeric array containing 'shape' and 'rate' parameters for the Gamma distribution. Both shape and rate should be > 0.")
  
  hyperpar = c(shape = 1, 1) 
  expect_identical(check_hyperpar_gamma(hyperpar), "hyperpar should be a named numeric array containing 'shape' and 'rate' parameters for the Gamma distribution. Both shape and rate should be > 0.")
  
  hyperpar = c(1, 1) 
  expect_identical(check_hyperpar_gamma(hyperpar), "hyperpar should be a named numeric array containing 'shape' and 'rate' parameters for the Gamma distribution. Both shape and rate should be > 0.")
  
  hyperpar = c(shape = 1, rate = -1)
  expect_identical(check_hyperpar_gamma(hyperpar), "hyperpar should be a named numeric array containing 'shape' and 'rate' parameters for the Gamma distribution. Both shape and rate should be > 0.")
  
})



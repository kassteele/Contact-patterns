context("final")

#if (getOption("priorSens_test_computefinal")){
  
test_that("'gaussian case' works as expected", {
  
  mm0<-70; ll0<-0.5; nn<-4; kk<-1; mh<-76
  
  n <-nn
  const<-rep(1,n)
  kappa<-kk
  sd.kappa<-sqrt(1/kappa)
  mu.dat<-mh
  set.seed(100)
  y = rnorm(n,mu.dat,sd.kappa)
  y<-y-mean(y)+mu.dat
  ex1.data<-data.frame(y,const)
  
  # Running INLA
  #--------------
  
  ex1.m0l0 = inla(y ~ const-1, data = ex1.data,
                  control.fixed=list(mean=list(const=mm0),prec=list(const=ll0)), 
                  control.family = list(hyper = list(prec = list(initial=log(kk),fixed=TRUE))),
                  num.threads=1)
  
  inla_marginal_posterior <- ex1.m0l0$"marginals.fixed"$"const" # marginal posterior density at (mm0,ll0)
  
  # Computing sensibility
  #-----------------------
  
  grid_epsilon <- getOption("priorSens_grid_epsilon") # epsilon used in the grid exploration around prior values.
  number_axis_points <- 200  # How many points for a given parameter in the grid exploration, 
  
  gaussian_grid_polar = compute_grid_polar(number_axis_points = 2*number_axis_points, 
                                           log_prior_density = "gaussian",
                                           hyperpar = c(mean = mm0, prec = ll0))
  
  inla_sens_gaussian = compute_sensitivity(grid = gaussian_grid_polar,
                                           inla_marginal_posterior = inla_marginal_posterior, 
                                           integration_limits = TRUE)
  
  max_inla_sens_gaussian = max(inla_sens_gaussian$sensitivity[,"sensitivity"])
  
  # Analytical computation
  analytical_sensitivity = compute_analytical_sensitivity_gaussian(grid_values = gaussian_grid_polar$cartesian, 
                                                                   prior_mean = mm0,
                                                                   prior_precision = ll0,
                                                                   sample_size = n,
                                                                   kappa = kappa, 
                                                                   sample_mean = mu.dat,
                                                                   grid_epsilon = grid_epsilon)
  
  max_analytical_sensitivity = max(analytical_sensitivity[,3])
  
  
  expect_equal(max_inla_sens_gaussian, max_analytical_sensitivity, tolerance = 1e-2)
  
  sensitivity_obj = prior_sensitivity(log_prior_density = "gaussian",
                                      hyperpar = c(mean = mm0, prec = ll0),
                                      inla_marginal_posterior = inla_marginal_posterior, 
                                      number_axis_points = 2*number_axis_points,
                                      grid_epsilon = grid_epsilon,                                      
                                      integration_limits = TRUE)
  
  expect_identical(sensitivity_obj$polar_grid, gaussian_grid_polar)
  expect_identical(sensitivity_obj$sensitivity, inla_sens_gaussian$sensitivity)
  expect_equal(max(sensitivity_obj$sensitivity[,"sensitivity"]), max_analytical_sensitivity, tolerance = 1e-2)

  # minimal call
  sensitivity_obj = prior_sensitivity(log_prior_density = "gaussian",
                                      hyperpar = c(mean = mm0, prec = ll0),
                                      inla_marginal_posterior = inla_marginal_posterior)
  
  expect_equal(max(sensitivity_obj$sensitivity[,"sensitivity"]), max_analytical_sensitivity, tolerance = 1e-2)

  # Testing user defined case
    
  log_f_normal <- function(x, hyperpar){
    dnorm(x, mean = hyperpar[1], sd = 1/sqrt(exp(hyperpar[2])), log = TRUE)
  }
  
  gaussian_grid_polar_user = compute_grid_polar(number_axis_points = 2*number_axis_points,
                                                log_prior_density = log_f_normal,
                                                hyperpar = c(mm0, log(ll0)),
                                                mean_prior_density = mm0,
                                                grid_epsilon = grid_epsilon)
  
  inla_sens_gaussian_user = compute_sensitivity(grid = gaussian_grid_polar_user,
                                                inla_marginal_posterior = inla_marginal_posterior)
  
  max_inla_sens_gaussian_user = max(inla_sens_gaussian_user$sensitivity[,"sensitivity"])
  expect_equal(max_inla_sens_gaussian_user, max_analytical_sensitivity, tolerance = 1e-2)
  
  inla_sens_gaussian_user2 = prior_sensitivity(number_axis_points = 2*number_axis_points,
                                               log_prior_density = log_f_normal,
                                               hyperpar = c(mm0, log(ll0)),
                                               mean_prior_density = mm0,
                                               inla_marginal_posterior = inla_marginal_posterior,
                                               grid_epsilon = grid_epsilon)
  
  max_inla_sens_gaussian_user = max(inla_sens_gaussian_user2$sensitivity[,"sensitivity"])
  expect_equal(max_inla_sens_gaussian_user, max_analytical_sensitivity, tolerance = 1e-2)
  
})

test_that("'gamma case' works as expected", {
  
  ## Parameters for the model in Ex2
  aa0<-0.5
  bb0<-0.01
  n <- 10
  m <- 2
  sigma2.ML<-0.035
  ntsigma2.ML<-n*sigma2.ML
  set.seed(100)
  y <- rnorm(n, mean = m, sd = sqrt(bb0 / aa0))
  y <- m + sqrt(sigma2.ML) * (y - mean(y)) / (sqrt((n - 1) / n) * sd(y))
  stopifnot(all.equal(ntsigma2.ML,
                      sum((y - m)^2)))
  ex2.data <- data.frame(y, const=1)
  ex2.a0b0 <- inla(y ~ -1 + const,
                   data = ex2.data,
                   control.fixed=list(mean=list(const=m),prec=list(const=1e12)),
                   control.family=list(hyper =list(prec =list(prior="loggamma",param=c(aa0, bb0)))),
                   control.compute=list(hyperpar=TRUE),
                   control.predictor = list(initial = 15),
                   num.threads=1)
  ex2.a0b0.hyp <- inla.hyperpar(ex2.a0b0, dz = 0.25, diff.logdens = 15)
  
  ex2.post.log.a0b0<-ex2.a0b0.hyp$"internal.marginals.hyperpar"$`Log precision for the Gaussian observations`
  exact_posterior <- cbind(ex2.post.log.a0b0[,1], dgamma(exp(ex2.post.log.a0b0[,1]), shape = aa0 + n/2, rate = bb0 + 0.5 * ntsigma2.ML) * exp(ex2.post.log.a0b0[,1]))
  
  # Computing sensibility
  #-----------------------
  
  number_axis_points = 100
  grid_epsilon = getOption("priorSens_grid_epsilon")
  
  gamma_grid_polar = compute_grid_polar(number_axis_points = 2*number_axis_points, 
                                        log_prior_density = "gamma",
                                        hyperpar = c(shape = aa0, rate = bb0),
                                        grid_epsilon = grid_epsilon)
  
  inla_sens_gamma = compute_sensitivity(grid = gamma_grid_polar,
                                        inla_marginal_posterior = ex2.post.log.a0b0, 
                                        integration_limits = TRUE)
  
  max_inla_sens_gamma = max(inla_sens_gamma$sensitivity[,"sensitivity"])
  
  # Analytical computation
  analytical_sensitivity = compute_analytical_sensitivity_gamma(grid_values = gamma_grid_polar$cartesian, 
                                                                prior_shape = aa0,
                                                                prior_rate = bb0,
                                                                sample_size = n,
                                                                sum_squares = sigma2.ML*n,
                                                                grid_epsilon = grid_epsilon)
  
  max_analytical_sensitivity = max(analytical_sensitivity[,3])  
  
  expect_equal(max_inla_sens_gamma, max_analytical_sensitivity, tolerance = 1e-3)
  
  # minimal call
  sensitivity_obj = prior_sensitivity(log_prior_density = "gamma",
                                      hyperpar = c(shape = aa0, rate = bb0),
                                      inla_marginal_posterior = ex2.post.log.a0b0)
  
  expect_equal(max(sensitivity_obj$sensitivity[,"sensitivity"]), max_analytical_sensitivity, tolerance = 1e-3)
  
  # Testing user defined gamma case
  
  log_loggamma <- function(x, hyperpar){
    dgamma(exp(x), shape = exp(hyperpar[1]), rate = exp(hyperpar[2]), log = TRUE) + x
  }
  
  gamma_user_grid_polar1 = compute_grid_polar(number_axis_points = 2*number_axis_points, 
                                              log_prior_density = log_loggamma,
                                              grid_epsilon = grid_epsilon,
                                              hyperpar = c(log(aa0), log(bb0)))
  
  inla_sens_gamma_user = compute_sensitivity(grid = gamma_user_grid_polar1,
                                             inla_marginal_posterior = ex2.post.log.a0b0)
  
  max_inla_sens_gamma_user = max(inla_sens_gamma_user$sensitivity[,"sensitivity"])
  
  expect_equal(max_inla_sens_gamma_user, max_analytical_sensitivity, tolerance = 1e-3)
  
  inla_sens_gamma_user2 = prior_sensitivity(number_axis_points = 2*number_axis_points,
                                                 log_prior_density = log_loggamma,
                                                 hyperpar = c(log(aa0), log(bb0)),
                                                 inla_marginal_posterior = ex2.post.log.a0b0)
  
  max_inla_sens_gamma_user = max(inla_sens_gamma_user2$sensitivity[,"sensitivity"])
  
  expect_equal(max_inla_sens_gamma_user, max_analytical_sensitivity, tolerance = 1e-3)
  
})

test_that("'rw1 case' works as expected", {
  
  set.seed(2524)
  
  # simulate data
  n0 <- 100
  x0 <- rep(NA, n0)
  prec.rw1 <- 3
  prec.y <- 0.5
  x0[1] <- 0
  for(i in 2:n0){
    x0[i] <- x0[i-1] + rnorm(1, mean=0, sd=sqrt(1/prec.rw1))
  }
  y <- x0 + rnorm(n0, mean=0, sd=sqrt(1/prec.y))
  
  # run inla
  data0 <- data.frame(y=y, x=1:n0)
  model <- y~f(x,model="rw1", hyper=list(prec=list(prior="loggamma",param=c(1, 0.005))), constr=F) -1
  result <- inla(model, family="gaussian", data=data0, 
                 control.family=list(hyper=list(prec=list(initial=log(prec.y), fixed=T))),
                 control.inla=list(strategy="laplace"),
                 control.predictor=list(compute=TRUE), 
                 verbose=F,
                 control.compute=list(hyperpar=TRUE),
                 num.threads=1)
  ex3.m <- inla.hyperpar(result)
  # INLA log marginal posterior needed for sensitivity computations
  ex3.m.log.prec.a0b0<-ex3.m$internal.marginals.hyperpar$"Log precision for x"
  
  # compute grid
  grid_epsilon <- 0.05
  number_axis_points = 200
  
  ex3_gamma_grid_polar = compute_grid_polar(number_axis_points = number_axis_points, 
                                            log_prior_density = "gamma",
                                            hyperpar = c(shape = 1, rate = 0.005),
                                            grid_epsilon = grid_epsilon)
  
  # compute sensitivity
  ex3_sensitivity = compute_sensitivity(grid = ex3_gamma_grid_polar,
                                        inla_marginal_posterior = ex3.m.log.prec.a0b0, 
                                        integration_limits = TRUE)
  
  max_ex3_sensitivity <- max(ex3_sensitivity$sensitivity[,"sensitivity"])   # polar sensitivity = 0.8334672                                               
  
  expect_equal(max_ex3_sensitivity, 0.8334672, tolerance = 1e-2)
  
  
})

#}
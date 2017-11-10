old.compute_hellinger_distance <- function # Computation of the Hellinger distance between two inla posteriors
(inla_posterior_1, 
 ### inla posterior marginal 1
 inla_posterior_2,
 ### inla posterior marginal 2
 integration_limits = FALSE,
 ### Set TRUE if you want to use 10% rule for integration limits, 
 ### it uses Inf by default.
 renormalize = TRUE
 ### If TRUE (default), both 'inla_posterior_1' and 'inla_posterior_2'
 ### will renormalized to guarantee they integrate to 1.
){

  if (renormalize){
    
    range_post_1 <- range(inla_posterior_1[,1])
    length_range_post_1 = range_post_1[2] - range_post_1[1]  
    
    cte1 = integrate(f=inla.dmarginal, 
                     range_post_1[1] - 0.1*length_range_post_1, 
                     range_post_1[2] + 0.1*length_range_post_1,  
                     marginal = inla_posterior_1)$value
    
    inla_posterior_1[,2] = inla_posterior_1[,2]/cte1
    
    range_post_2 <- range(inla_posterior_2[,1])
    length_range_post_2 = range_post_2[2] - range_post_2[1]  
    
    cte2 = integrate(f=inla.dmarginal, 
                     range_post_2[1] - 0.1*length_range_post_2, 
                     range_post_2[2] + 0.1*length_range_post_2,  
                     marginal = inla_posterior_2)$value
    
    inla_posterior_2[,2] = inla_posterior_2[,2]/cte2
    ##details<<
    ## Before computing the hellinger distance, we renormalize both posteriors just
    ## to make sure they integrate to 1, this turns out to be inefficient if we
    ## already know they are properly normalized. In this case set 'renormalize = FALSE'
    
  }

  bc_integrand <- cbind(inla_posterior_1[,1], sqrt(inla_posterior_1[,2]*inla_posterior_2[,2]))
  
  range_bc_integrand <- range(bc_integrand[,1])
  length_range = range_bc_integrand[2] - range_bc_integrand[1]
  
  integration_limit_inf = - Inf
  integration_limit_sup = Inf
  
  if (integration_limits){
    integration_limit_inf = range_bc_integrand[1] - 0.1*length_range
    integration_limit_sup = range_bc_integrand[2] + 0.1*length_range
  }
  
  bc_coef <- integrate(inla.dmarginal, 
                       integration_limit_inf, 
                       integration_limit_sup, 
                       marginal = bc_integrand)$value
  
  bc_coef = min(bc_coef, 1)
  
  return(sqrt(1-bc_coef))
  ##value<<
  ## return the Hellinger distance.
  
}

#----------------------------------------------------------------------------------------------------------------

# Functions useful to compute the Gaussian grid and sensitivity
#---------------------------------------------------------------

compute_precision_range_gaussian <- function # Compute range of precision in the Gaussian case
### This function compute the range of the precision parameter that has roots
### in the grid exploration used to compute prior sensitivity.
(prior_precision,
 ### prior precision used in the analysis.
 grid_epsilon
 ### epsilon used in the definition of the contour
 ){

# Q_MR_5: I have replaced "eps" in the code by grid_epsilon in order
# to make the code more consistent. Do you agree?
  
  temp1 = sqrt(1-(1-grid_epsilon^2)^4)
  temp2 = ((1-grid_epsilon^2)^4)
  
  precision_inf <- prior_precision*(1-temp1)^2/temp2
  precision_sup <- prior_precision*(1+temp1)^2/temp2
  
  return(c(precision_inf, precision_sup))
  ##value<<
  ## Return the endpoints of the precision range.
  
}

compute_mean_roots_gaussian <- function # Compute mean roots for a given new precision in the Gaussian case.
### This function computes analytically the two mean roots for each new precision value
### in the sensibility grid exploration of the Gaussian case.
(new_precision, 
 ### a new precision value in the contour.
 prior_mean, 
 ### prior mean used in the analysis.
 prior_precision, 
 ### prior precision used in the analysis.
 grid_epsilon
 ### epsilon used in the definition of the contour
 ){

# Q_MR_6: I have replaced "eps" in the code by grid_epsilon in order
# to make the code more consistent. Do you agree?
  
  temp1 = prior_precision+new_precision
  temp2 = prior_precision*new_precision
  temp3 = sqrt(abs(-4*(temp1)*log(sqrt((temp1)/(2*sqrt(temp2)))*(1-grid_epsilon^2))/(temp2)))
  ##details<<
  ## Absolute value under sqrt is taken to prevent numerically extremely 
  ## small but NEGATIVE values under it.
    
  mean_inf <- prior_mean - temp3
  mean_sup <- prior_mean + temp3

# Q_MR_7: Do you agree that return(c(mean_inf, mean_sup, 
# rep(new_precision,2))) is used? This was we will have mean values 
# first and precision values afterwards. For me it is a more natural
# way of parametrisation of the normal distribution.
  
  return(c(rep(new_precision,2), mean_inf, mean_sup))
  ##value<<
  ## For each new_precision value, it returns a vector with 4 elements,
  ## first the precisions and then the means.
  
}

old.compute_grid_gaussian <- function # Grid computation for a Gaussian prior
(number_axis_points,
 ### Number of precision points used in the grid.
 prior_mean, 
 ### prior mean used in the analysis.
 prior_precision, 
 ### prior precision used in the analysis.
 grid_epsilon
 ### epsilon used in the definition of the contour.
){
  
  precision_range <- compute_precision_range_gaussian(prior_precision, grid_epsilon)

  axis.points <- seq(precision_range[1], precision_range[2], length.out = number_axis_points)

  grid_roots <- compute_mean_roots_gaussian(axis.points, prior_mean, prior_precision, grid_epsilon)
  grid_roots <- matrix(grid_roots, length(grid_roots)/2, 2)
  colnames(grid_roots) <- c("precision", "mean")

  return(grid_roots)
  ##value<<
  ## returns a matrix with two columns, the first with precision values and 
  ## the second with mean values

}

compute_corrected_posteriors_gaussian <- function # Compute corrected posteriors in the Gaussian case.
### This function computes corrected posteriors for values of the epsilon-grid in the Gaussian case
(new_mean, 
 ### new prior mean
 new_precision, 
 ### new precision mean
 inla_marginal_posterior, 
 ### posterior marginal returned by INLA.
 prior_mean, 
 ### prior mean used in the INLA analysis.
 prior_precision,
 ### prior precision used in the INLA analysis
 integration_limits = FALSE
 ### Set TRUE if you want to use 10% rule for integration limits, 
 ### it uses Inf by default.
){
  
  new_prior <- dnorm(inla_marginal_posterior[,1], mean = new_mean, sd = 1/sqrt(new_precision))
  base_prior <- dnorm(inla_marginal_posterior[,1], mean = prior_mean, sd = 1/sqrt(prior_precision))
  
  corrected_posterior_marginal <- inla_marginal_posterior[,2] * new_prior/base_prior
  
  unnormalized_corrected_posterior_marginal <- cbind(inla_marginal_posterior[,1], corrected_posterior_marginal)
  
  range_marginal_posterior <- range(unnormalized_corrected_posterior_marginal[,1])
  length_range = range_marginal_posterior[2] - range_marginal_posterior[1]
  
  integration_limit_inf = - Inf
  integration_limit_sup = Inf
  
  if (integration_limits){
    integration_limit_inf = range_marginal_posterior[1] - 0.1*length_range
    integration_limit_sup = range_marginal_posterior[2] + 0.1*length_range
  }
  
  normalization_constant <- integrate(inla.dmarginal, 
                                      integration_limit_inf, 
                                      integration_limit_sup, 
                                      marginal = unnormalized_corrected_posterior_marginal)$value
  
  corrected_posterior_marginal = cbind(unnormalized_corrected_posterior_marginal[,1],
                                       unnormalized_corrected_posterior_marginal[,2]/normalization_constant)
  
  return(corrected_posterior_marginal)
  ##value<<
  ## Corrected posterior marginals using INLA format.
  
}

old.compute_sensitivity_gaussian <- function # Compute sensitivity measures for every point in the Grid
(grid_values, 
 ### Grid values computed with 'old.compute_grid_gaussian' function.
 prior_mean, 
 ### prior mean used in the analysis.
 prior_precision, 
 ### prior precision used in the analysis.
 inla_marginal_posterior, 
 ### posterior marginal returned by INLA.
 grid_epsilon, 
 ### epsilon used in the definition of the contour
 integration_limits = FALSE
 ### Set TRUE if you want to use 10% rule for integration limits, 
 ### it uses Inf by default.
){
  
  grid_size = nrow(grid_values)
  sensibility_vector <- rep(NA, grid_size)
  
  for (i in 1:grid_size){
    
    new_mean <- grid_values[i, "mean"]
    new_precision <- grid_values[i, "precision"]
    
    corrected_posterior_marginal = compute_corrected_posteriors_gaussian(new_mean, new_precision, 
                                                                         inla_marginal_posterior, prior_mean, prior_precision,
                                                                         integration_limits = integration_limits)
    
    sensibility_vector[i] <- old.compute_hellinger_distance(inla_posterior_1 = inla_marginal_posterior, 
                                                        inla_posterior_2 = corrected_posterior_marginal,
                                                        integration_limits)/grid_epsilon
    
  }
  
  sensibility_matrix <- cbind(grid_values, sensibility = sensibility_vector)
  
  return(sensibility_matrix)
  
}

#----------------------------------------------------------------------------------------------------------------

# Functions useful to compute the Gamma grid and sensitivity
#-------------------------------------------------------------

old.compute_hellinger_gammas <- function # Hellinger distance between two Gammas
(prior_shape, 
 ### shape parameter of the base prior
 prior_rate, 
 ### rate parameter of the base prior
 new_shape, 
 ### shape parameter of the new prior
 new_rate, 
 ### rate parameter of the new prior
 analytic = TRUE
 ### Use analytical computation? Default is TRUE.
 ){
  
  if (!analytic){ # The option to compute the HD numerically is very unstable
                  # when used inside 'uniroot.all. Need more testing, for now
                  # use only the option 'analytic = TRUE'.
# Q_MR_9: I am very sorry to hear that. Is this difficulty possibly
# connected to the extremely wide interval = c(0, sup_interval_shape)
# for uniroot.all? 
    
    g <- function(x, prior_shape, prior_rate, new_shape, new_rate){
      exp(0.5 *(dgamma(x, shape = prior_shape, rate = prior_rate, log = TRUE) + 
                  dgamma(x, shape = new_shape, rate = new_rate, log = TRUE))) 
    }
# Q_MR_10: Maybe we could use the log-Gamma representation leading to
    #    g <- function(x, prior_shape, prior_rate, new_shape, new_rate){
    #      exp(x + 0.5 *(dgamma(exp(x), shape = prior_shape, rate = prior_rate, log = TRUE) + 
    #                  dgamma(exp(x), shape = new_shape, rate = new_rate, log = TRUE))) 
    #    }
# Q_MR_11: If you agree with Q_MR_10 the lower integration range
# should be set to lower = -Inf
    bc_coef = integrate(f = g, lower = 0, upper = Inf, 
                        prior_shape = prior_shape, 
                        prior_rate = prior_rate, 
                        new_shape = new_shape, 
                        new_rate = new_rate)$value
    bc_coef = min(bc_coef, 1)
    
    hellinger_distance = sqrt(1 - bc_coef)
    return(hellinger_distance)
    
  } 

  a0 = prior_shape; a1 = new_shape; b0 = prior_rate; b1 = new_rate
  
  hellinger_distance = sqrt(1 - gamma((a0 + a1)/2)*sqrt((b0^a0 * b1^a1)/(gamma(a0) * gamma(a1) * ((b0 + b1)/2)^(a0+a1))))
  return(hellinger_distance)
  
}

compute_roots_shape_given_rate <- function # Compute shape parameter values given a fixed rate parameter.
### This function returns shape parameters on the epsilon grid for a given rate parameter
(new_rate, 
 ### The rate parameter for which you want to find shape parameter values on the epsilon grid.
 prior_shape, 
 ### shape parameter of the base prior.
 prior_rate, 
 ### rate parameter of the base prior.
 grid_epsilon, 
 ### epsilon used on the epsilon grid computations.
 max_num_iter = 10,
 ### If the roots are not contained on the initial interval, we double the 
 ### range of the interval in each iteration.
 analytic = TRUE
 ){
# Q_MR_12: Is it necessary for the function to be "hidden" by starting
# the name by a dot? What is the purpose of this point in front of the
# function name?
  .f_fixed_rate <- function(x, prior_shape, prior_rate, fixed_rate, grid_epsilon){
    
    old.compute_hellinger_gammas(prior_shape = prior_shape, 
                             prior_rate = prior_rate, 
                             new_shape = x, 
                             new_rate = fixed_rate,
                             analytic = analytic) - grid_epsilon
    
  }
  f_fixed_rate = Vectorize(.f_fixed_rate, vectorize.args="x") # Find the roots of this function
  
  sup_interval_shape = 2*prior_shape
# Q_MR_13: Deals with all interval settings below: Is interval = c(0,
# sup_interval_shape) too wide for uniroot.all to give precise
# results?
  shape_range = try(uniroot.all(f = f_fixed_rate, interval = c(0, sup_interval_shape), 
                            prior_shape = prior_shape, 
                            prior_rate = prior_rate, 
                            fixed_rate = new_rate, 
                            grid_epsilon = grid_epsilon))
  
  if (inherits(shape_range, "try-error")) shape_range = rep(NA, 2) 
  
  num_iter = 1
  
  while ((length(shape_range) < 2) & (num_iter < max_num_iter)){ 
    
    sup_interval_shape = 2*sup_interval_shape # Double the interval
    
    shape_range = try(uniroot.all(f = f_fixed_rate, interval = c(0, sup_interval_shape), 
                              prior_shape = prior_shape, 
                              prior_rate = prior_rate, 
                              fixed_rate = new_rate, 
                              grid_epsilon = grid_epsilon))
    
    if (inherits(shape_range, "try-error")) shape_range = rep(NA, 2)
    
    num_iter = num_iter + 1
    
  }
  
  if (length(shape_range) != 2) shape_range = rep(NA, 2) 
  # Sometimes it finds more than 2 roots due to numerical instability on the computation of
  # the Hellinger distance (when it is numerically computed instead of analytically computed)
  return(shape_range)
  
}

compute_roots_rate_given_shape <- function # Compute rate parameter values given a fixed shape parameter.
### This function returns rate parameters on the epsilon grid for a given shape parameter
(new_shape, 
 ### The shape parameter for which you want to find rate parameter values on the epsilon grid. 
 prior_shape, 
 ### shape parameter of the base prior.
 prior_rate, 
 ### rate parameter of the base prior.
 grid_epsilon, 
 ### epsilon used on the epsilon grid computations.
 max_num_iter = 10,
 ### If the roots are not contained on the initial interval, we double the 
 ### range of the interval in each iteration.
 analytic = TRUE
 ){
  
  .f_fixed_shape <- function(x, prior_shape, prior_rate, fixed_shape, grid_epsilon){
    
    old.compute_hellinger_gammas(prior_shape = prior_shape, 
                             prior_rate = prior_rate, 
                             new_shape = fixed_shape, 
                             new_rate = x,
                             analytic = analytic) - grid_epsilon
    
  }
  f_fixed_shape = Vectorize(.f_fixed_shape, vectorize.args="x") # Funtion to find the root of
  
  sup_interval_rate = 2*prior_rate
  
  rate_range = uniroot.all(f = f_fixed_shape, interval = c(0, sup_interval_rate), 
                            prior_shape = prior_shape, 
                            prior_rate = prior_rate, 
                            fixed_shape = new_shape, 
                            grid_epsilon = grid_epsilon)
  
  if (inherits(rate_range, "try-error")) rate_range = rep(NA, 2)
  
  num_iter = 1
  
  while ((length(rate_range) < 2) & (num_iter < max_num_iter)){
    
    sup_interval_rate = 2*sup_interval_rate # Double the interval
    
    rate_range = uniroot.all(f = f_fixed_shape, interval = c(0, sup_interval_rate), 
                              prior_shape = prior_shape, 
                              prior_rate = prior_rate, 
                              fixed_shape = new_shape, 
                              grid_epsilon = grid_epsilon)
    
    if (inherits(rate_range, "try-error")) rate_range = rep(NA, 2)
    
    num_iter = num_iter + 1
    
  }
  
  if (length(rate_range) != 2) rate_range = rep(NA, 2)
  # Sometimes it finds more than 2 roots due to numerical instability on the computation of
  # the Hellinger distance (when it is numerically computed instead of analytically computed)
  
  return(rate_range)
  
}

old.compute_grid_gamma <- function # Compute the epsilon grid for the Gamma prior
(number_axis_points, 
 ### Number of points used in the grid.
 prior_shape, 
 ### shape parameter of the base prior.
 prior_rate, 
 ### rate parameter of the base prior.
 grid_epsilon, 
 ### epsilon used in the definition of the contour.
 max_num_iter = 10,
 ### If the roots are not contained on the initial interval, we double the 
 ### range of the interval in each iteration.
 analytic = TRUE
 ){
  
  shape_range =  suppressWarnings(compute_roots_shape_given_rate(new_rate = prior_rate,
                                                                 prior_shape, 
                                                                 prior_rate, 
                                                                 grid_epsilon, 
                                                                 max_num_iter = max_num_iter,
                                                                 analytic = analytic))
  
  axis_shape_points <- seq(shape_range[1], shape_range[2], length.out = number_axis_points)

  
  gamma_grid1 = cbind(rep(axis_shape_points, 2),
                     as.numeric(t(suppressWarnings(sapply(axis_shape_points, FUN = compute_roots_rate_given_shape, 
                                                          prior_shape = prior_shape, 
                                                          prior_rate = prior_rate, 
                                                          grid_epsilon = grid_epsilon, 
                                                          max_num_iter = max_num_iter,
                                                          analytic = analytic)))))
  colnames(gamma_grid1) <- c("shape", "rate")
  
  axis_rate_points <- unique(gamma_grid1[,"rate"])
  
  gamma_grid2 = cbind(as.numeric(t(suppressWarnings(sapply(axis_rate_points, FUN = compute_roots_shape_given_rate, 
                                                          prior_shape = prior_shape, 
                                                          prior_rate = prior_rate, 
                                                          grid_epsilon = grid_epsilon, 
                                                          max_num_iter = max_num_iter,
                                                           analytic = analytic)))), rep(axis_rate_points, 2))
  colnames(gamma_grid2) <- c("shape", "rate")
  
  gamma_grid2 = gamma_grid2[!is.na(gamma_grid2[,"shape"]), ]
  
  gamma_grid = rbind(gamma_grid1, gamma_grid2)
  return(gamma_grid)
  
}

compute_corrected_posteriors_gamma <- function # Compute corrected posteriors in the Gamma case.
### This function computes corrected posteriors for values of the epsilon-grid in the Gamma case
(new_shape, 
 ### new shape parameter
 new_rate, 
 ### new rate parameter
 inla_marginal_posterior, 
 ### posterior marginal returned by INLA.
# Q_MR_14: Should we mention that only internal marginal posteriors
# (in log-scale) can be used here? Do we need to check if the 
# inla_marginal_posterior is in the correct scale?
 
 prior_shape, 
 ### shape parameter of the base prior.
 prior_rate,
 ### rate parameter of the base prior.
 integration_limits = FALSE
 ### Set TRUE if you want to use 10% rule for integration limits, 
 ### it uses Inf by default.
){
  
  new_prior <- dgamma(exp(inla_marginal_posterior[,1]), shape = new_shape, rate = new_rate) * exp(inla_marginal_posterior[,1])
  base_prior <- dgamma(exp(inla_marginal_posterior[,1]), shape = prior_shape, rate = prior_rate) * exp(inla_marginal_posterior[,1])
  
  corrected_posterior_marginal <- inla_marginal_posterior[,2] * new_prior/base_prior
  
  unnormalized_corrected_posterior_marginal <- cbind(inla_marginal_posterior[,1], corrected_posterior_marginal)
  
  range_marginal_posterior <- range(unnormalized_corrected_posterior_marginal[,1])
  length_range = range_marginal_posterior[2] - range_marginal_posterior[1]
  
  integration_limit_inf = - Inf
  integration_limit_sup = Inf
  
  if (integration_limits){
    integration_limit_inf = range_marginal_posterior[1] - 0.1*length_range
    integration_limit_sup = range_marginal_posterior[2] + 0.1*length_range
  }
  
  normalization_constant <- integrate(inla.dmarginal, 
                                      integration_limit_inf, 
                                      integration_limit_sup, 
                                      marginal = unnormalized_corrected_posterior_marginal)$value
  
  corrected_posterior_marginal = cbind(unnormalized_corrected_posterior_marginal[,1],
                                       unnormalized_corrected_posterior_marginal[,2]/normalization_constant)
  
  return(corrected_posterior_marginal)
  ##value<<
  ## Corrected posterior marginals using INLA format.
  
}

old.compute_sensitivity_gamma <- function # Compute sensitivity measures for every point in the Grid
(grid_values, 
 ### Grid value scomputed with 'old.compute_grid_gamma' function.
 prior_shape, 
 ### shape parameter of the base prior.
 prior_rate, 
 ### rate parameter of the base prior.
 inla_marginal_posterior, 
 ### posterior marginal returned by INLA.
 grid_epsilon, 
 ### epsilon used in the definition of the contour
 integration_limits = FALSE
 ### Set TRUE if you want to use 10% rule for integration limits, 
 ### it uses Inf by default.
){
  
  grid_size = nrow(grid_values)
  sensibility_vector <- rep(NA, grid_size)
  
  for (i in 1:grid_size){
    
    new_shape <- grid_values[i, "shape"]
    new_rate <- grid_values[i, "rate"]
    
    corrected_posterior_marginal = compute_corrected_posteriors_gamma(new_shape, new_rate, 
                                                                      inla_marginal_posterior, 
                                                                      prior_shape, prior_rate,
                                                                      integration_limits = integration_limits)
    
        
    sensibility_vector[i] <- old.compute_hellinger_distance(inla_posterior_1 = inla_marginal_posterior, 
                                                        inla_posterior_2 = corrected_posterior_marginal,
                                                        integration_limits)/grid_epsilon
    
  }
  
  sensibility_matrix <- cbind(grid_values, sensibility = sensibility_vector)
  
  return(sensibility_matrix)
  
}

old.compute_analytical_sensitivity_gamma_ex2 <- function # Compute sensitivity measures for every point in the Grid
(grid_values, 
 ### Grid values computed with 'old.compute_grid_gamma' function.
 prior_shape, 
 ### shape parameter of the base prior.
 prior_rate, 
 ### rate parameter of the base prior.
 sample_size,
 ### sample size of the experiment.
 sum_squares,
 ### sum((y - mean)^2)
 grid_epsilon, 
 ### epsilon used in the definition of the contour
 integration_limits = FALSE
 ### Set TRUE if you want to use 10% rule for integration limits, 
 ### it uses Inf by default.
){
  
  grid_size = nrow(grid_values)
  sensibility_vector <- rep(NA, grid_size)
  
  for (i in 1:grid_size){
    
    new_shape <- grid_values[i, "shape"]
    new_rate <- grid_values[i, "rate"]

    sensibility_vector[i] <- old.compute_hellinger_gammas(prior_shape = prior_shape + sample_size/2, 
                                                      prior_rate = prior_rate + 0.5 * sum_squares, 
                                                      new_shape = new_shape + sample_size/2, 
                                                      new_rate = new_rate + 0.5 * sum_squares, 
                                                      analytic = TRUE)/grid_epsilon
    
  }
  
  sensibility_matrix <- cbind(grid_values, sensibility = sensibility_vector)
  
  return(sensibility_matrix)
  
}

#----------------------------------------------------------------------------------------------------------------



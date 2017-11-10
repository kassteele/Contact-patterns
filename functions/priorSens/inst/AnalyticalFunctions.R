https://www.dropbox.com/home/Sensitivity%20project/Sensitivity20130725/package_additional_function...eps05R64.R


Ex1:

old.compute_hellinger_gaussians with analytic=TRUE

old.compute_analytical_sensitivity_gaussian_ex1_polar


Ex2:

old.compute_hellinger_gammas with analytic=TRUE

old.compute_analytical_sensitivity_gamma_ex2



#----------------------------------------------------------------------------------------------------------------

# Functions useful to compute the Gaussian sensitivity in Ex1 analytically
#---------------------------------------------------------------


old.compute_hellinger_gaussians <- function # Hellinger distance between two Gaussians
(prior_mean, 
 ### mean parameter of the base prior
 prior_precision, 
 ### precision parameter of the base prior
 new_mean, 
 ### mean parameter of the new prior
 new_precision, 
 ### precision parameter of the new prior
 analytic = TRUE
 ### Use analytical computation? Default is TRUE.
){
  
  if (!analytic){ # MR: Not tested yet.
    
    
    g <- function(x, prior_mean, prior_precision, new_mean, new_precision){
      exp(0.5 *(dnorm(x, mean = prior_mean, sd = 1/sqrt(prior_precision), log = TRUE) + 
                  dnorm(x, mean = new_mean, sd = 1/sqrt(new_precision), log = TRUE))) 
    }
    
    bc_coef = integrate(f = g, lower = -Inf, upper = Inf, 
                        prior_mean = prior_mean, 
                        prior_precision = prior_precision, 
                        new_mean = new_mean, 
                        new_precision = new_precision)$value
    
    hellinger_distance = sqrt(1 - bc_coef)
    return(hellinger_distance)
    
  } 
  
  m0 = prior_mean; m1 = new_mean; l0 = prior_precision; l1 = new_precision
  
  hellinger_distance = sqrt(1 - sqrt((2/sqrt(l1*l0))/(1/l1 + 1/l0))*exp(-(m1 - m0)^2/(4*(1/l1 + 1/l0))))
  return(hellinger_distance)
  
}



old.compute_analytical_sensitivity_gaussian_ex1_polar <- function # Compute sensitivity measures for every point in the Grid
(grid_values, 
 ### Grid values computed with 'compute_grid_gaussian' function.
 prior_mean, 
 ### prior mean used in the analysis.
 prior_precision, 
 ### prior precison used in the analysis.
 sample_size,
 ### sample size of the experiment.
 kappa, 
 ### kappa
 sample_mean,
 ### sample mean
 grid_epsilon, 
 ### epsilon used in the definition of the contour
 integration_limits = FALSE
 ### Set TRUE if you want to use 10% rule for integration limits, 
 ### it uses Inf by default.
){
  
  grid_size = nrow(grid_values)
  sensibility_vector <- rep(NA, grid_size)
  
  for (i in 1:grid_size){
    
    new_precision <- grid_values[i, "prec"]
    new_mean <- grid_values[i, "mean"]
    
    sensibility_vector[i] <- old.compute_hellinger_gaussians(prior_mean = (sample_size*kappa*sample_mean+prior_precision*prior_mean)/(sample_size*kappa+prior_precision), 
                                                             prior_precision = sample_size*kappa+prior_precision, 
                                                             new_mean = (sample_size*kappa*sample_mean+new_precision*new_mean)/(sample_size*kappa+new_precision), 
                                                             new_precision = sample_size*kappa+new_precision, 
                                                             analytic = TRUE)/grid_epsilon
    
  }
  
  sensibility_matrix <- cbind(grid_values, sensitivity = sensibility_vector)
  
  return(sensibility_matrix)
  
}


#----------------------------------------------------------------------------------------------------------------

# Functions useful to compute the Gaussian sensitivity in Ex2 analytically
#---------------------------------------------------------------


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
  
  sensibility_matrix <- cbind(grid_values, sensitivity = sensibility_vector)
  
  return(sensibility_matrix)
  
}
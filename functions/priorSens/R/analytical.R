# Functions only useful for testing purposes
#--------------------------------------------

compute_hellinger_gaussians_analytic <- function # Hellinger distance between two Gaussians
(prior_mean, 
 ### mean parameter of the base prior
 prior_precision, 
 ### precision parameter of the base prior
 new_mean, 
 ### mean parameter of the new prior
 new_precision
 ### precision parameter of the new prior
){
  
  m0 = prior_mean; m1 = new_mean; l0 = prior_precision; l1 = new_precision
  
  hellinger_distance = sqrt(1 - sqrt((2/sqrt(l1*l0))/(1/l1 + 1/l0))*exp(-(m1 - m0)^2/(4*(1/l1 + 1/l0))))
  return(hellinger_distance)
  
}

#-----

compute_analytical_sensitivity_gaussian <- function # Compute analytical sensitivity measures for every point in the grid
(grid_values, 
 ### cartesian grid values
 prior_mean, 
 ### prior mean used in the analysis.
 prior_precision, 
 ### prior precison used in the analysis.
 sample_size,
 ### sample size of the experiment.
 kappa, 
 ### fixed observation precision
 sample_mean,
 ### sample mean
 grid_epsilon
 ### epsilon used in the definition of the contour
){
  
  grid_size = nrow(grid_values)
  sensibility_vector <- rep(NA, grid_size)
  
  for (i in 1:grid_size){
    
    new_precision <- grid_values[i, "prec"]
    new_mean <- grid_values[i, "mean"]
    
    base_posterior_mean = (sample_size*kappa*sample_mean+prior_precision*prior_mean)/(sample_size*kappa+prior_precision)
    base_posterior_prec = sample_size*kappa+prior_precision
    
    new_posterior_mean = (sample_size*kappa*sample_mean+new_precision*new_mean)/(sample_size*kappa+new_precision)
    new_posterior_prec = sample_size*kappa+new_precision
    
    sensibility_vector[i] <- compute_hellinger_gaussians_analytic(prior_mean = base_posterior_mean, 
                                                                  prior_precision = base_posterior_prec, 
                                                                  new_mean = new_posterior_mean, 
                                                                  new_precision = new_posterior_prec)/grid_epsilon
    
  }
  
  sensibility_matrix <- cbind(grid_values, sensitivity = sensibility_vector)
  
  return(sensibility_matrix)
  
}

#--------------------------------------------------------------------------------------------------------------

compute_hellinger_gammas_analytic <- function # Hellinger distance between two Gammas
(prior_shape, 
 ### shape parameter of the base prior
 prior_rate, 
 ### rate parameter of the base prior
 new_shape, 
 ### shape parameter of the new prior
 new_rate 
 ### rate parameter of the new prior
){
  
  a0 = prior_shape; a1 = new_shape; b0 = prior_rate; b1 = new_rate
  
  hellinger_distance = sqrt(1 - gamma((a0 + a1)/2)*sqrt((b0^a0 * b1^a1)/(gamma(a0) * gamma(a1) * ((b0 + b1)/2)^(a0+a1))))
  return(hellinger_distance)
  
}

#-----

compute_analytical_sensitivity_gamma <- function # Compute analytical sensitivity measures for every point in the grid
(grid_values, 
 ### cartesian grid
 prior_shape, 
 ### shape parameter of the base prior.
 prior_rate, 
 ### rate parameter of the base prior.
 sample_size,
 ### sample size of the experiment.
 sum_squares,
 ### sum((y - mean)^2)
 grid_epsilon
 ### epsilon used in the definition of the contour
){
  
  grid_size = nrow(grid_values)
  sensibility_vector <- rep(NA, grid_size)
  
  for (i in 1:grid_size){
    
    new_shape <- grid_values[i, "shape"]
    new_rate <- grid_values[i, "rate"]
    
    sensibility_vector[i] <- compute_hellinger_gammas_analytic(prior_shape = prior_shape + sample_size/2, 
                                                               prior_rate = prior_rate + 0.5 * sum_squares, 
                                                               new_shape = new_shape + sample_size/2, 
                                                               new_rate = new_rate + 0.5 * sum_squares)/grid_epsilon
    
  }
  
  sensibility_matrix <- cbind(grid_values, sensitivity = sensibility_vector)
  
  return(sensibility_matrix)
  
}
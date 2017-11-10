
# Functions to compute corrected posteriors
#-------------------------------------------

compute_ratio_priors_gaussian <- function # Compute the ratio 'new_prior/base_prior' for Gaussian case.
(hyperpar_new, 
 ### Hyperparameter of the new prior. Needs to be a named numerical array with 'mean' and 
 ### 'prec' parameters.
 hyperpar_old, 
 ### Hyperparameter of the old/base prior. Needs to be a named numerical array with 'mean' and 
 ### 'prec' parameters.
 inla_internal_marginal
 ### Needs to be the inla posterior marginal on the theta scale where theta is assigned a 
 ### Gaussian prior.
){
  
  new_mean = as.numeric(hyperpar_new[1]); new_prec = as.numeric(hyperpar_new[2])
  prior_mean = as.numeric(hyperpar_old[1]); prior_prec = as.numeric(hyperpar_old[2])
  
  new_prior <- dnorm(inla_internal_marginal[,1], mean = new_mean, sd = 1/sqrt(new_prec))
  base_prior <- dnorm(inla_internal_marginal[,1], mean = prior_mean, sd = 1/sqrt(prior_prec))
  
  return(new_prior/base_prior)
  
}

compute_ratio_priors_gamma <- function # Compute the ratio 'new_prior/base_prior' for gamma case.
(hyperpar_new, 
 ### Hyperparameter of the new prior. Needs to be a named numerical array with 'shape' 
 ### and 'rate' parameters.
 hyperpar_old, 
 ### Hyperparameter of the old/base prior. Needs to be a named numerical array with 'shape' 
 ### and 'rate' parameters.
 inla_internal_marginal
 ### Needs to be the inla posterior marginal on the theta scale where theta is assigned a 
 ### Gamma prior.
){
  
  new_shape = as.numeric(hyperpar_new[1]); new_rate = as.numeric(hyperpar_new[2])
  prior_shape = as.numeric(hyperpar_old[1]); prior_rate = as.numeric(hyperpar_old[2])
  
  log_new_prior <- dgamma(exp(inla_internal_marginal[,1]), shape = new_shape, rate = new_rate, log = TRUE) + inla_internal_marginal[,1]
  log_base_prior <- dgamma(exp(inla_internal_marginal[,1]), shape = prior_shape, rate = prior_rate, log = TRUE) + inla_internal_marginal[,1]
  
  return(exp(log_new_prior - log_base_prior))
  
}

compute_ratio_priors_custom <- function # Compute the ratio 'new_prior/base_prior' for custom case.
(hyperpar_new, 
 ### Hyperparameter of the new prior. Needs to be a numerical array.
 hyperpar_old, 
 ### Hyperparameter of the old/base prior. Needs to be a numerical array.
 inla_internal_marginal,
 ### Needs to be the inla posterior marginal on the theta scale where theta is assigned a 
 ### Gamma prior.
 log_prior_density
 ### For user defined prior distributions, it should be a 
 ### function that takes 'x' defined on the 
 ### real line and a numeric vector of length 2 called 'hypepar' with the hyperparameters of 
 ### 'log_prior_density' also defined on the real line as arguments.
 ### log_prior_density(x, hypepar) must return the log density value evaluated at 'x'.
){
  
  log_new_prior <- log_prior_density(inla_internal_marginal[,1], hyperpar = hyperpar_new)
  log_base_prior <- log_prior_density(inla_internal_marginal[,1], hyperpar = hyperpar_old)
  
  return(exp(log_new_prior - log_base_prior))
  
}

compute_corrected_posteriors <- function # Compute corrected posteriors from INLA marginals.
### This function computes corrected posteriors for hyperparameters values of the epsilon-grid given
### an INLA posterior marginal in the appropriate representation.
(prior_type, 
 ### A character vector indication what is the family of the prior. Currently, available for 
 ### "gaussian", "gamma" and "custom".  
 hyperpar_new,
 ### hyperparameters of the new prior, see Details. 
 hyperpar_old,
 ### hyperparameters of the old/base prior, see Details. 
 inla_marginal_posterior, 
 ### Needs to be the inla posterior marginal on the theta scale where theta is assigned a 
 ### 'prior_type'.
 method = getOption("priorSens_correc_post_integration_method"),
 ### integration method used to compute the normalization cte of the corrected posteriors.
 ### Currently, the following methods are available: 
 ### 'integrate', 'simpson13rule'. Default is 'simpson13rule'.
 integration_limits = TRUE,
 ### Set TRUE if you want to use 10% rule for integration limits, otherwise uses Inf. Default is TRUE.
 log_prior_density = NULL
 ### For user defined prior distributions, it should be a 
 ### function that takes 'x' defined on the 
 ### real line and a numeric vector of length 2 called 'hypepar' with the hyperparameters of 
 ### 'log_prior_density' also defined on the real line as arguments.
 ### log_prior_density(x, hypepar) must return the log density value evaluated at 'x'.
){
  
  ##details<<
  ## Currently, this function is implemented only for 'prior_type' equal to 'gaussian' and 'gamma'. 
  msg_prior_type = "compute_corrected_posteriors: Currently 'prior_type' needs to be 'gaussian', 'gamma' or 'custom'."
  if (is.null(prior_type)){
    stop(msg_prior_type)
  } else if (!(prior_type %in% c("gaussian", "gamma", "custom"))){
    stop(msg_prior_type)
  }
  
  ##details<<
  ## If prior_type="gaussian", this includes a named numerical array with 'mean' and 'prec' 
  ## parameters. If prior_type="gamma", this includes a named numerical array with
  ## 'shape' and 'rate' parameters.
  if (prior_type == "gaussian"){
    ratio_prior = compute_ratio_priors_gaussian(hyperpar_new, hyperpar_old, inla_internal_marginal = inla_marginal_posterior)
  } else if (prior_type == "gamma"){
    ratio_prior = compute_ratio_priors_gamma(hyperpar_new, hyperpar_old, inla_internal_marginal = inla_marginal_posterior)
  } else if (prior_type == "custom"){
    ratio_prior = compute_ratio_priors_custom(hyperpar_new, hyperpar_old, inla_internal_marginal = inla_marginal_posterior,
                                              log_prior_density = log_prior_density)
  } else {
    stop(msg_prior_type)
  }
  
  corrected_posterior_marginal <- inla_marginal_posterior[,2] * ratio_prior
  
  unnormalized_corrected_posterior_marginal <- cbind(inla_marginal_posterior[,1], corrected_posterior_marginal)
  
  if (method == "integrate"){
    
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
                                        rel.tol = getOption("priorSens_integrate_rel.tol"),
                                        abs.tol = getOption("priorSens_integrate_abs.tol"),
                                        subdivisions = getOption("priorSens_integrate_subdivisions"),
                                        marginal = unnormalized_corrected_posterior_marginal)$value
  } else if (method == "simpson13rule"){
    
    
    normalization_constant <- simp13rule_matrix(density_matrix = unnormalized_corrected_posterior_marginal)
    
  }
  
  corrected_posterior_marginal = cbind(unnormalized_corrected_posterior_marginal[,1],
                                       unnormalized_corrected_posterior_marginal[,2]/normalization_constant)
  
  return(corrected_posterior_marginal)
  ##value<<
  ## Corrected posterior marginals using INLA format.
  
}

#----------------------------------------------------------------------------------------------------------------

compute_sensitivity <- function # Compute sensitivity measures for every point in the Grid
(grid, 
 ### Grid values computed with 'compute_grid_polar' function.
 inla_marginal_posterior, 
 ### posterior marginal returned by INLA.
 integration_correc_post_method = getOption("priorSens_correc_post_integration_method"),
 ### integration method used to compute the normalization cte of the corrected posteriors.
 ### Currently, the following methods are available: 
 ### 'integrate', 'simpson13rule'. Default is 'simpson13rule'.
 integration_inla_method = getOption("priorSens_inla_HD_integration_method"),
 ### integration method used to within the 'compute_hellinger_distance_inla'
 ### Currently, the following methods are available: 
 ### 'integrate', 'simpson13rule'. Default is 'simpson13rule'.
 integration_limits = TRUE
 ### Set TRUE if you want to use 10% rule for integration limits, otherwise uses Inf. Default is TRUE.
){
  
  if (!inherits(grid, "grid"))
    stop("'grid' argument needs to be an object with class 'grid'.'")
  
  if (!is.matrix(inla_marginal_posterior)){
    stop("'inla_marginal_posterior' is supposed to be a matrix.")
  } else if (ncol(inla_marginal_posterior) != 2)
    stop("'inla_marginal_posterior' is supposed to be a matrix with 2 columns.")
  
  grid_epsilon = grid$grid_epsilon
  
  grid_cartesian = grid$cartesian
  grid_size = nrow(grid_cartesian)
  
  grid_polar = grid$polar
  
  base_hyperpar = hyperpar_new = grid$hyperpar
  #names_hyperpar = names(base_hyperpar)
  
  prior_type = grid$prior_type
  log_prior_density = grid$log_prior_density
  
  # Normalizing inla posterior marginal
  cte = simp13rule_matrix(density_matrix = inla_marginal_posterior)
  inla_marginal_posterior[,2] = inla_marginal_posterior[,2]/cte
  
  sensitivity_vector <- rep(NA, grid_size)
  
  for (i in 1:grid_size){
    
    hyperpar_new = as.numeric(grid_cartesian[i,])
    
    corrected_posterior_marginal = compute_corrected_posteriors(prior_type = prior_type, 
                                                                hyperpar_new = hyperpar_new, 
                                                                hyperpar_old = base_hyperpar,
                                                                inla_marginal_posterior = inla_marginal_posterior, 
                                                                method = integration_correc_post_method, 
                                                                integration_limits = integration_limits,
                                                                log_prior_density = log_prior_density)
    
    sensitivity_vector[i] <- compute_hellinger_distance_inla(inla_posterior_1 = inla_marginal_posterior, 
                                                             inla_posterior_2 = corrected_posterior_marginal,
                                                             renormalize = FALSE,
                                                             method = integration_inla_method,
                                                             integration_limits)/grid_epsilon
    
  }
  
  sensitivity_matrix <- cbind(grid_polar, grid_cartesian, sensitivity = sensitivity_vector)
  
  result <- list(polar_grid = grid,
                 sensitivity = sensitivity_matrix)
  
  return(result)
  
}

#----------------------------------------------------------------------------------------------------------------

prior_sensitivity <- function # Compute grid and prior sensitivity values.
### Computes both the grid and the sensitvity values of a given
### prior distribution.
(log_prior_density,
 ### For Gaussian priors use 'gaussian', for gamma priors use 'gamma'.
 ### For user defined prior distributions, it should be a 
 ### function that takes 'x' defined on the 
 ### real line and a numeric vector of length 2 called 'hypepar' with the hyperparameters of 
 ### 'log_prior_density' also defined on the real line as arguments.
 ### log_prior_density(x, hypepar) must return the log density value evaluated at 'x'.
 hyperpar,
 ### For the Gaussian case needs to be a named array with 'mean' and 'prec' elements.
 ### For the Gamma case needs to be a named array with 'shape' and 'rate' elements.
 ### For the user defined case, needs to be a numeric vector that is consistent with 
 ### 'log_prior_density'.
 inla_marginal_posterior, 
 ### posterior marginal returned by INLA.
 number_axis_points = getOption("priorSens_grid_number_axis_points"),
 ### Number of points used in the grid.
 grid_epsilon = getOption("priorSens_grid_epsilon"),
 ### epsilon used in the definition of the contour.
 method = getOption("priorSens_grid_root_method"),
 ### Methods used to find the radius in the polar coordinate system. 
 ### Currently, the following methods are available: 
 ### 'nlminb', 'uniroot.all', 'uniroot'. Default is 'uniroot'.
 parallel = getOption("priorSens_grid_parallel"),
 ### Use parallel computing through the multicore package.
 coord_explr_method = getOption("priorSens_grid_coord_explr_method"),
 ### This is advanced option and sets the method used to find the root
 ### under coordinate exploration mode, which is more delicate and 
 ### therefore needs to use a more robust method.
 coord_explr_method_options = getOption("priorSens_grid_coord_explr_method_options"),
 ### A list with options to be passed to 'method_coord_exploration' routine.
 integration_grid_method = getOption("priorSens_grid_HD_integration_method"),
 ### Integration method in case we use numerical integration to compute HD.
 domain = getOption("priorSens_grid_HD_integration_domain"),
 ### Necessary only when computing the HD numerically. The end-points of the hellinger 
 ### distance integral.
 number_points = getOption("priorSens_grid_HD_integration_number_points"),
 ### Number of points used in the integration routine when simpson13rule is used.
 integration_correc_post_method = getOption("priorSens_correc_post_integration_method"),
 ### integration method used to compute the normalization cte of the corrected posteriors.
 ### Currently, the following methods are available: 
 ### 'integrate', 'simpson13rule'. Default is 'simpson13rule'.
 integration_inla_method = getOption("priorSens_inla_HD_integration_method"),
 ### integration method used to within the 'compute_hellinger_distance_inla'
 ### Currently, the following methods are available: 
 ### 'integrate', 'simpson13rule'. Default is 'simpson13rule'.
 integration_limits = TRUE,
 ### Set TRUE if you want to use 10% rule for integration limits, otherwise uses Inf. Default is TRUE.
 ...
 ### Extra parameters for 'compute_roots_r_given_thetas' function.
){
  
  polar_grid = compute_grid_polar(number_axis_points = number_axis_points, 
                                  log_prior_density = log_prior_density,
                                  hyperpar = hyperpar,
                                  grid_epsilon = grid_epsilon,
                                  method = method,
                                  coord_explr_method = coord_explr_method,
                                  coord_explr_method_options = coord_explr_method_options,
                                  integration_method = integration_grid_method,
                                  domain = domain,
                                  number_points = number_points,
                                  parallel = parallel,
                                  ...)
  
  sensitivity = compute_sensitivity(grid = polar_grid, 
                                    inla_marginal_posterior = inla_marginal_posterior,
                                    integration_correc_post_method = integration_correc_post_method,
                                    integration_inla_method = integration_inla_method,
                                    integration_limits = integration_limits)
  
  return(sensitivity)
  
}
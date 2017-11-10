# Functions useful to compute grid
#------------------------------------

.coord_correc_func <- function # Specify which coordinate correction to use.
(theta, 
 ### theta parameter of the polar coordinate, within [-pi, pi].
 coordinate_correction
 ### NULL if no coordinate correction is to be applied. Otherwise, a list
 ### with the following elements: 'xneg', 'xpos', 'yneg', 'ypos' indicating
 ### the coordinate corrections that will be applied dependent on the value of theta, 
 ### see Details.
){
  
  ##details<<
  ## The theta space is divided into 4 quadrants: if theta belongs to 
  ## [-pi, -pi/2] then we are in the negative side of the x-xis and of the y-axis.
  ## if theta belongs to [-pi/2, 0] we are on the positive side of the x-axis and on
  ## the negative side of the y axis, and so on.
  
  if (is.null(coordinate_correction)){ # in case both axis have a very different magnitude
    coord_correc_x = 1; coord_correc_y = 1; 
  } else if (theta <= -pi/2){
    coord_correc_x = coordinate_correction$xneg; coord_correc_y = coordinate_correction$yneg;      
  } else if (theta <= 0){
    coord_correc_x = coordinate_correction$xpos; coord_correc_y = coordinate_correction$yneg;      
  } else if (theta <= pi/2){
    coord_correc_x = coordinate_correction$xpos; coord_correc_y = coordinate_correction$ypos;      
  } else {
    coord_correc_x = coordinate_correction$xneg; coord_correc_y = coordinate_correction$ypos;      
  }
  
  return(c(coord_correc_x, coord_correc_y))
  
}
coord_correc_func <- cmpfun(.coord_correc_func)

#------------------------------------------------------------------------------

.map_hyperpar_gaussian <- function # Map from polar coordinate to cartesian coordinate for Gaussian prior.
(hyperpar,
 ### The initial hyperparameters of the Gaussian prior. Needs to be a numeric 
 ### vector where the first element is the mean and the second element is the precision.
 r, 
 ### radius of the polar coordinate.
 theta, 
 ### The degree of the polar coordinate.
 coord_correc_x, 
 ### x-axis correction, see function 'coord_correc_func' for details.
 coord_correc_y
 ### y-axis correction, see function 'coord_correc_func' for details.
){
  
  new_mean = hyperpar[1] + (r * cos(theta))*coord_correc_x
  log_new_prec =  log(hyperpar[2]) + (r * sin(theta))*coord_correc_y
  
  result = cbind(mean = new_mean, prec = exp(log_new_prec))
  return(result)
  
}
map_hyperpar_gaussian <- cmpfun(.map_hyperpar_gaussian)

objective_to_compute_roots_r_given_theta_gaussian <- function # Return HD - epsilon between two Gaussian priors
(x, 
 ### x = log(r), where r is the radius.
 hyperpar, 
 ### should be a numeric vector where the first element contain the mean and the second element contain the precision.
 fixed_theta, 
 ### a given theta within [-pi, pi]
 grid_epsilon, 
 ### epsilon used on the epsilon grid computations.
 coordinate_correction,
 ### this will specify adjustments necessary to go from polar to cartesian coordinates
 ### in order to cover the whole space efficiently. See function 'coord_correc_func' for details.
 log_prior_density = NULL,
 ### log prior density function. Not necessary for the 'gaussian' case.
 ...
){
  
  ########################
  #count <<- count + 1
  ########################
  
  coord_correc = coord_correc_func(theta = fixed_theta, coordinate_correction)
  
  new_hyperpar = map_hyperpar_gaussian(hyperpar = hyperpar, 
                                       r = exp(x), 
                                       theta = fixed_theta, 
                                       coord_correc_x = coord_correc[1], 
                                       coord_correc_y = coord_correc[2])
  
  new_mean = as.numeric(new_hyperpar[,"mean"])
  new_prec =  as.numeric(new_hyperpar[,"prec"])
  ### Since x belongs to the real line, r = exp(x) > 0 as it should be.
  
  compute_hellinger(log_dist1 = "gaussian", hyperpar1 = cbind(mean = hyperpar[1], prec = hyperpar[2]),
                    log_dist2 = "gaussian", hyperpar2 = cbind(mean = new_mean, prec = new_prec),
                    domain = NULL) - grid_epsilon
  
}

#---------------------------------------------------------------------------------------------------------

.map_hyperpar_gamma <- function # Map from polar coordinate to cartesian coordinate for Gamma prior.
(hyperpar,
 ### The initial hyperparameters of the Gamma prior. Needs to be a numeric 
 ### vector where the first element is the shape and the second element is the rate.
 r, 
 ### radius of the polar coordinate.
 theta, 
 ### The degree of the polar coordinate.
 coord_correc_x, 
 ### x-axis correction, see function 'coord_correc_func' for details.
 coord_correc_y
 ### y-axis correction, see function 'coord_correc_func' for details.
){
  
  log_new_shape = log(hyperpar[1]) + (r * cos(theta))*coord_correc_x
  log_new_rate =  log(hyperpar[2]) + (r * sin(theta))*coord_correc_y
  
  result = cbind(shape = exp(log_new_shape), rate = exp(log_new_rate))
  return(result)
  
}
map_hyperpar_gamma <- cmpfun(.map_hyperpar_gamma)

objective_to_compute_roots_r_given_theta_gamma <- function # Return HD - epsilon between two Gamma priors
(x, 
 ### x = log(r), where r is the radius.
 hyperpar, 
 ### should be a named array with 'shape' and 'rate' elements.
 fixed_theta, 
 ### a given theta within [-pi, pi]
 grid_epsilon, 
 ### epsilon used on the epsilon grid computations.
 coordinate_correction,
 ### this will specify adjustments necessary to go from polar to cartesian coordinates
 ### in order to cover the whole space efficiently. See function 'coord_correc_func' for details.
 log_prior_density = NULL,
 ### log prior density function. Not necessary for the 'gamma' case.
 ...
){
  
  coord_correc = coord_correc_func(theta = fixed_theta, coordinate_correction)
  
  new_hyperpar = map_hyperpar_gamma(hyperpar = hyperpar, 
                                    r = exp(x), 
                                    theta = fixed_theta, 
                                    coord_correc_x = coord_correc[1], 
                                    coord_correc_y = coord_correc[2])
  
  new_shape = new_hyperpar[,"shape"]
  new_rate =  new_hyperpar[,"rate"]
  ### Since x belongs to the real line, r = exp(x) > 0 as it should be.
  
  compute_hellinger(log_dist1 = "gamma", hyperpar1 = cbind(shape = hyperpar[1], rate = hyperpar[2]),
                    log_dist2 = "gamma", hyperpar2 = cbind(shape = new_shape, rate = new_rate),
                    domain = NULL) - grid_epsilon
  
}

#---------------------------------------------------------------------------------------------------------

.map_hyperpar_custom <- function # Map from polar coordinate to cartesian coordinate for an user-defined prior.
(hyperpar,
 ### The initial hyperparameters, defined on the real line, of the user-defined prior. Needs to be a numeric 
 ### vector of length 2.
 r, 
 ### radius of the polar coordinate.
 theta, 
 ### The degree of the polar coordinate.
 coord_correc_x, 
 ### x-axis correction, see function 'coord_correc_func' for details.
 coord_correc_y
 ### y-axis correction, see function 'coord_correc_func' for details.
){
  
  new_hyper1 = hyperpar[1] + (r * cos(theta))*coord_correc_x
  new_hyper2 = hyperpar[2] + (r * sin(theta))*coord_correc_y
  
  result = cbind(new_hyper1, new_hyper2)
  return(result)
  
}
map_hyperpar_custom = cmpfun(.map_hyperpar_custom)

objective_to_compute_roots_r_given_theta_custom <- function # Return HD - epsilon between two user-defined priors
(x, 
 ### x = log(r), where r is the radius.
 hyperpar, 
 ### should be a numeric vector, where the elements can take value on the real line.
 fixed_theta, 
 ### a given theta within [-pi, pi]
 grid_epsilon, 
 ### epsilon used on the epsilon grid computations.
 coordinate_correction,
 ### this will specify adjustments necessary to go from polar to cartesian coordinates
 ### in order to cover the whole space efficiently. See function 'coord_correc_func' for details.
 log_prior_density = NULL,
 ### log prior density function. Should be a function that takes 'x' defined on the 
 ### real line and a numeric vector of length 2 called 'hypepar' with the hyperparameters of 
 ### 'log_prior_density' also defined on the real line as arguments.
 ### log_prior_density(x, hypepar) must return the log density value evaluated at 'x'.
 integration_method = "integrate",
 ### Integration method in case we use numerical integration.
 domain = c(-Inf, Inf),
 ### Necessary only when computing the HD numerically. The end-points of the hellinger 
 ### distance integral.
 number_points = 100
 ### Number of points used in the integration routine when simpson13rule is used.
){
  
  if (grid_epsilon <= 0)
    stop("'grid_epsilon' needs to be greater than zero.")
  
  coord_correc = coord_correc_func(theta = fixed_theta, coordinate_correction)
  
  coord_correc_x = coord_correc[1]
  coord_correc_y = coord_correc[2]
  
  prior_1 = as.numeric(hyperpar[1])
  prior_2 = as.numeric(hyperpar[2])
  
  new_hyperpar = map_hyperpar_custom(hyperpar = hyperpar, 
                                    r = exp(x), 
                                    theta = fixed_theta, 
                                    coord_correc_x = coord_correc_x, 
                                    coord_correc_y = coord_correc_y)
  
  new_prior_1 = as.numeric(new_hyperpar[1])
  new_prior_2 = as.numeric(new_hyperpar[2])
  ### Since x belongs to the real line, r = exp(x) > 0 as it should be.
  
  compute_hellinger(log_dist1 = log_prior_density, hyperpar1 = c(prior_1, prior_2),
                    log_dist2 = log_prior_density, hyperpar2 = c(new_prior_1, new_prior_2),
                    integration_method = integration_method,
                    domain = domain,
                    number_points = number_points) - grid_epsilon
    
}

#------------------------------------------------------------------------------

compute_roots_r_given_thetas <- function # Compute radial parameter values given a fixed theta.
### This function returns radial parameters on the epsilon grid for a given theta parameter. 
(new_theta, 
 ### The theta parameter for which you want to find radial parameter values on the epsilon grid.
 ### Needs to be between [-pi, pi].
 objective_function,
 ### objective function used to compute the roots. Should take the form 
 ### f(x, hyperpar, fixed_theta, grid_epsilon, coordinate_correction, log_prior_density) and
 ### return (HD - grid_epsilon). There are also some pre-defined cases, like
 ### objective_to_compute_roots_r_given_theta_gaussian for Gaussian, 
 ### objective_to_compute_roots_r_given_theta_gamma for gamma and
 ### objective_to_compute_roots_r_given_theta_custom for user-defined prior functions.
 hyperpar,
 ### For the Gaussian case needs to be a named array with 'mean' and 'prec' elements.
 ### For the Gamma case needs to be a named array with 'shape' and 'rate' elements.
 ### For a user defined function should be consistent to what the defined function requires, 
 ### but needs to be a numeric vector of length 2 defined on the real line.
 grid_epsilon,
 ### epsilon used on the epsilon grid computations.
 coordinate_correction = NULL,
 ### this will specify adjustments necessary to go from polar to cartesian coordinates
 ### in order to cover the whole space efficiently. See function 'coord_correc_func' for details.
 log_prior_density = NULL,
 ### log prior density function. For user defined prior distributions, it should be a 
 ### function that takes 'x' defined on the 
 ### real line and a numeric vector of length 2 called 'hypepar' with the hyperparameters of 
 ### 'log_prior_density' also defined on the real line as arguments.
 ### log_prior_density(x, hypepar) must return the log density value evaluated at 'x'.
 ### This argument should be NULL for particular cases already implemented, like 'gaussian'
 ### and 'gamma'. 
 method = getOption("priorSens_grid_root_method"),
 ### Methods used to find the radius of the polar coordinate system. 
 ### Currently, the following methods are available: 'nlminb', 'uniroot.all', 'uniroot'. 
 ### Default is 'uniroot'.
 search_interval = c(-3, 2),
 ### Only if method equals 'uniroot.all' or 'uniroot'. A vector containing the end-points 
 ### of the interval to be searched for the root of r in log scale. If no root in found, 
 ### try to increase this interval. Remember that we have applied some kind of orthogonalization
 ### in the polar coordinate system, so that r is expected to be around 1, so a 
 ### search_interval = c(-3, 3) implies a search interval (0.049, 20.085) for r, which is likely
 ### to contain the root.
 integration_method = "integrate",
 ### Integration method in case we use numerical integration.
 domain = c(-Inf, Inf),
 ### Necessary only when computing the HD numerically. The end-points of the hellinger 
 ### distance integral.
 number_points = 100,
 ### Number of points used in the integration routine when simpson13rule is used.
 ...
 ### Extra arguments for the 'method' used.
){
  
  if ((new_theta < -pi) | (new_theta > pi)) stop("new_theta needs to be between [-pi, pi]")
  
  if (grid_epsilon <= 0)
    stop("'grid_epsilon' needs to be greater than zero.")
  
  if ((!is.null(coordinate_correction)) & 
        (!all(c("xneg", "xpos", "yneg", "ypos") %in% names(coordinate_correction))))
    stop("'coordinate_correction' needs to be a list with the following ellements 'xneg', 'xpos', 'yneg', 'ypos'")
  
  f_fixed_theta <- match.fun(objective_function)
  f_fixed_theta_square <- function(x, ...) f_fixed_theta(x, ...)^2
  
  if (method == "nlminb"){
    
    new_r = suppressWarnings(try(nlminb(start = 0, # Because I expect r = exp(x) to be close to 1
                                        objective = f_fixed_theta_square, 
                                        hyperpar = hyperpar, 
                                        fixed_theta = new_theta, 
                                        grid_epsilon = grid_epsilon, 
                                        coordinate_correction = coordinate_correction,
                                        log_prior_density = log_prior_density,
                                        integration_method = integration_method,
                                        domain = domain,
                                        number_points = number_points,
                                        ...)$par, 
                                 silent = TRUE))
    
    
  } else if (method == "uniroot.all"){
    
    new_r = suppressWarnings(try(uniroot.all(f = f_fixed_theta, interval = search_interval, 
                                             hyperpar = hyperpar, fixed_theta = new_theta, 
                                             grid_epsilon = grid_epsilon, 
                                             coordinate_correction = coordinate_correction,
                                             log_prior_density = log_prior_density, 
                                             integration_method = integration_method,
                                             domain = domain,
                                             number_points = number_points,
                                             ...), 
                                 silent = TRUE))
    
  } else if (method == "uniroot"){
    
    new_r = suppressWarnings(try(uniroot(f = f_fixed_theta, interval = search_interval, 
                                         hyperpar = hyperpar, fixed_theta = new_theta, 
                                         grid_epsilon = grid_epsilon, 
                                         coordinate_correction = coordinate_correction,
                                         log_prior_density = log_prior_density, 
                                         integration_method = integration_method,
                                         domain = domain,
                                         number_points = number_points,
                                         ...)$root, 
                                 silent = TRUE))
    
    
  } else {stop("Invalid method.")}
  
  if (inherits(new_r, "try-error")){
    
    warning(new_r[1])
    new_r = NA
    
  } else if (length(new_r) != 1){
    
    warning("more than one solution found for radius, please check")
    new_r = NA
    
  }
  
  return(c(new_theta, exp(new_r)))
  
}

#------------------------------------------------------------------------------

compute_grid_polar <- function # Compute the epsilon grid using polar coordinates
(number_axis_points = getOption("priorSens_grid_number_axis_points"), 
 ### Number of points used in the grid.
 log_prior_density,
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
 mean_prior_density = NULL,
 ### In case 'log_prior_density' is a user-defined function,
 ### please provide the mean of the prior density of the parameter x 
 ### defined on the real line. 
 grid_epsilon = getOption("priorSens_grid_epsilon"),
 ### epsilon used in the definition of the contour.
 parallel = getOption("priorSens_grid_parallel"),
 ### Use parallel computing through the multicore package.
 method = getOption("priorSens_grid_root_method"),
 ### Methods used to find the radius in the polar coordinate system. 
 ### Currently, the following methods are available: 
 ### 'nlminb', 'uniroot.all', 'uniroot'. Default is 'uniroot'.
 coord_explr_method = getOption("priorSens_grid_coord_explr_method"),
 ### This is advanced option and sets the method used to find the root
 ### under coordinate exploration mode, which is more delicate and 
 ### therefore needs to use a more robust method.
 coord_explr_method_options = getOption("priorSens_grid_coord_explr_method_options"),
 ### A list with options to be passed to 'method_coord_exploration' routine.
 integration_method = getOption("priorSens_grid_HD_integration_method"),
 ### Integration method in case we use numerical integration to compute HD. For user-defined 
 ### cases the 'integration_method' will be set to "integrate"
 domain = getOption("priorSens_grid_HD_integration_domain"),
 ### Necessary only when computing the HD numerically. The end-points of the hellinger 
 ### distance integral.
 quantile_domain = getOption("priorSens_grid_HD_integration_domain_quantile"),
 ### In case 'domain = "compute_domain"', the integration limitis will be defined
 ### by the following quantiles: quantile_domain and (1 - quantile_domain)
 number_points = getOption("priorSens_grid_HD_integration_number_points"),
 ### Number of points used in the integration routine when simpson13rule is used.
 ...
 ### Extra parameters for 'compute_roots_r_given_thetas' function.
){
  
  if (number_axis_points <= 0) stop("'number_axis_points' needs to be > 0.")
  
  thetas_coord = c(-pi, -pi/2, 0, pi/2) # To get possible necessary correction factors 
  # for the coordinates.
  
  if (is.function(log_prior_density)){
    objective = match.fun(objective_to_compute_roots_r_given_theta_custom)
    prior_type = "custom"
    map_hyperpar = match.fun(map_hyperpar_custom)
    log_prior_density = match.fun(log_prior_density)
    if (any(domain == "compute_domain")){
      quantile1 = find_quantile(quantile = quantile_domain, log_dist = log_prior_density, 
                                hyperpar = hyperpar, mean_dist = mean_prior_density)
      quantile2 = find_quantile(quantile = 1 - quantile_domain, log_dist = log_prior_density, 
                                hyperpar = hyperpar, mean_dist = mean_prior_density)
      domain = c(quantile1, quantile2)
    }
    if (integration_method != "integrate")
      integration_method = "integrate"

    } else if (log_prior_density == "gaussian") {
    msg1 = check_hyperpar_gaussian(hyperpar = hyperpar)
    if (msg1 != "ok"){
      stop(paste("compute_grid_polar:", msg1))
    }
    hyperpar = c(as.numeric(hyperpar["mean"]), as.numeric(hyperpar["prec"]))
    objective = match.fun(objective_to_compute_roots_r_given_theta_gaussian)
    map_hyperpar = match.fun(map_hyperpar_gaussian)
    prior_type = "gaussian"
    log_prior_density <- NULL
  } else if (log_prior_density == "gamma") {
    msg1 = check_hyperpar_gamma(hyperpar = hyperpar)
    if (msg1 != "ok"){
      stop(paste("compute_grid_polar:", msg1))
    }
    objective = match.fun(objective_to_compute_roots_r_given_theta_gamma)
    map_hyperpar = match.fun(map_hyperpar_gamma)
    prior_type = "gamma"
    log_prior_density <- NULL
  } else {
    stop("Invalid 'log_prior_density'.")  
  } 
  
  grid_polar_coord = t(sapply(thetas_coord, 
                              FUN = compute_roots_r_given_thetas, 
                              objective_function = objective,
                              hyperpar = hyperpar,
                              grid_epsilon = grid_epsilon,
                              coordinate_correction = NULL,
                              log_prior_density = log_prior_density,
                              integration_method = integration_method,
                              domain = domain,
                              number_points = number_points,
                              method = coord_explr_method,
                              coord_explr_method_options))
  
  coordinate_correction = list(xneg = grid_polar_coord[which(grid_polar_coord[,1] == -pi), 2],
                               xpos = grid_polar_coord[which(grid_polar_coord[,1] == 0), 2],
                               yneg = grid_polar_coord[which(grid_polar_coord[,1] == -pi/2), 2],
                               ypos = grid_polar_coord[which(grid_polar_coord[,1] == pi/2), 2])
  
  thetas = seq(-pi, pi, length.out = number_axis_points)
  
  if (parallel){

    if (!require(multicore)) stop("You need multicore package to use 'parallel = TRUE' option.")
    
    grid_polar = matrix(unlist(mclapply(thetas, FUN = compute_roots_r_given_thetas, mc.silent=TRUE,
                                        objective_function = objective,
                                        hyperpar = hyperpar,
                                        grid_epsilon = grid_epsilon,
                                        coordinate_correction = coordinate_correction,
                                        log_prior_density = log_prior_density,
                                        integration_method = integration_method,
                                        domain = domain,
                                        number_points = number_points,
                                        method = method,
                                        ...)), number_axis_points, 2, byrow = T)
    
  } else {
    
    grid_polar = matrix(unlist(lapply(thetas, FUN = compute_roots_r_given_thetas, 
                                      objective_function = objective,
                                      hyperpar = hyperpar,
                                      grid_epsilon = grid_epsilon,
                                      coordinate_correction = coordinate_correction,
                                      log_prior_density = log_prior_density,
                                      integration_method = integration_method,
                                      domain = domain,
                                      number_points = number_points,
                                      method = method,
                                      ...)), number_axis_points, 2, byrow = T)
    
  }
  grid_polar_temp = grid_polar_coord
  grid_polar_temp[,2] = 1
  grid_polar = rbind(grid_polar, grid_polar_temp) # Add 4 points used in the coordinate exploration
                                                  # They have r = 1 by definition.
  
  grid_polar = grid_polar[order(grid_polar[,1]), ] # ordering
  colnames(grid_polar) <- c("theta", "r")
  
  if (any(is.na(grid_polar[, "r"]))){
    
    stop("NA found on the grid exploration, try different options for the optimization routine.")
    
  }
  
  coord_correc_thetas = t(sapply(grid_polar[,"theta"], 
                                 FUN=coord_correc_func, 
                                 coordinate_correction = coordinate_correction))
  
  grid_cartesian = map_hyperpar(hyperpar = hyperpar, 
                                r = grid_polar[,"r"], 
                                theta = grid_polar[,"theta"], 
                                coord_correc_x = unlist(coord_correc_thetas[,1]), 
                                coord_correc_y = unlist(coord_correc_thetas[,2]))
  
  result = list(polar = grid_polar, cartesian = grid_cartesian, hyperpar = hyperpar, 
                grid_epsilon = grid_epsilon, prior_type = prior_type, 
                log_prior_density = log_prior_density)
  class(result) <- "grid"
  return(result)
  
}


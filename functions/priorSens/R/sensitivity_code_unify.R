# Auxiliary functions
#----------------------

check_hyperpar_gaussian <- function(hyperpar){
  
  message = "ok"
  
  if (length(hyperpar) != 2){
    
    message = "hyperpar should be a named numeric array containing 'mean' and 'prec' parameters for the Gaussian distribution. 'prec' should be > 0."
    return(message)
    
  } else if ((is.null(names(hyperpar))) | (!all(names(hyperpar) == c("mean", "prec")))){
    
    message = "hyperpar should be a named numeric array containing 'mean' and 'prec' parameters for the Gaussian distribution. 'prec' should be > 0."
    return(message)
    
  } else if (hyperpar["prec"] <= 0){
    
    message = "hyperpar should be a named numeric array containing 'mean' and 'prec' parameters for the Gaussian distribution. 'prec' should be > 0."
    return(message)
    
  }
  
  return(message)
  
}

check_hyperpar_gamma <- function(hyperpar){

  message = "ok"
  
  if (length(hyperpar) != 2){
    
    message = "hyperpar should be a named numeric array containing 'shape' and 'rate' parameters for the Gamma distribution. Both shape and rate should be > 0."
    return(message)
    
  } else if ((is.null(names(hyperpar))) | (!all(names(hyperpar) == c("shape", "rate"))) | 
        (!all(hyperpar > 0))){
    
   message = "hyperpar should be a named numeric array containing 'shape' and 'rate' parameters for the Gamma distribution. Both shape and rate should be > 0."
   return(message)
   
  }
  
  return(message)

}

#---------------------------------------------------------------------------------------------------------------

# Integration routines
#----------------------

simp13rule_matrix <- function(density_matrix, apply_spline = TRUE){
  
  if (apply_spline)
    density_matrix = inla.smarginal(marginal=as.matrix(density_matrix), keep.type=TRUE)
  
  n = nrow(density_matrix)
  i.0 = c(1, n)
  i.4 = seq(2, n - 1, by = 2)
  i.2 = seq(3, n - 1, by = 2)
  
  integral = sum(sum(density_matrix[i.0, 2]) + 4*sum(density_matrix[i.4, 2]) + 2*sum(density_matrix[i.2, 2]))
  integral = integral*(diff(range(density_matrix[,1]))/(3*(n-1)))

  return(integral)

}

#---------------------------------------------------------------------------------------------------------------

# Functions to compute Hellinger distance
#------------------------------------------

compute_hellinger <- function # Hellinger distance between two distributions
### Computes the hellinger distance between two distributions numerically.
### The computation is made analitically in case both distributions are Gaussians
### or Gammas.
(log_dist1,
 ### should be a function that takes 'x' defined in the real line and 'hypepar' as arguments 
 ### and log_dist1(x, hyper) returns the log density value evaluated at 'x'. Write "gaussian" 
 ### or "gamma" in case both 'dist1' and 'dist2' are Gaussian or gamma distributions respectively.
 hyperpar1,
 ### extra parameter to be passed to 'log_dist1', as hyperparameters for example. 
 ### If "gaussian", this includes a named numerical array with 'mean' and 'prec' 
 ### parameters. If "gamma", this includes a named numerical array with
 ### 'shape' and 'rate' parameters.
 log_dist2,
 ### same as 'log_dist1'. 
 hyperpar2,
 ### same as 'hyperpar1'
 integration_method = "integrate",
 ### Integration method in case we use numerical integration.
 domain = c(-Inf, Inf),
 ### Necessary only when computing the HD numerically. The end-points of the hellinger 
 ### distance integral.
 number_points = 100
 ### Number of points used in the integration routine when simpson13rule is used.
){
  
  if (is.function(log_dist1) | is.function(log_dist2)){

    return(compute_hellinger_numerically_grid_pertubation(log_dist = log_dist1,
                                                          hyperpar1 = hyperpar1,
                                                          hyperpar2 = hyperpar2,
                                                          method = integration_method,
                                                          domain = domain,
                                                          number_points = number_points))
    
    # return(compute_hellinger_numerically(log_dist1, hyperpar1, log_dist2, hyperpar2, domain))      
    
  } else if ((log_dist1 == "gaussian") & (log_dist2 == "gaussian")){ # Gaussians case
    
    return(compute_hellinger_gaussians(hyperpar1, hyperpar2))
    
  } else if ((log_dist1 == "gamma") & (log_dist2 == "gamma")){ # Gammas case
    
    return(compute_hellinger_gammas(hyperpar1, hyperpar2))
    
  } else {
    
    stop("Invalid argument.")
    
  }
  
}

compute_hellinger_gaussians <- function # Compute Gaussian HD analitically
(hyperpar1, 
 ### Should be a named numeric array with 'mean' and 'prec' elements of 
 ### the first Gaussian distribution.
 hyperpar2
 ### Should be a named numeric array with 'mean' and 'prec' elements of 
 ### the second Gaussian distribution.
){
  
  #msg1 = check_hyperpar_gaussian(hyperpar1)
  #msg2 = check_hyperpar_gaussian(hyperpar2)
  #if ((msg1 != "ok") | (msg2 != "ok")){
  #  stop(paste("compute_hellinger_gaussians:", msg1))
  #}
  
  m1 = hyperpar1[,"mean"]; p1 = hyperpar1[,"prec"]
  m2 = hyperpar2[,"mean"]; p2 = hyperpar2[,"prec"]
  
  prodp = p1*p2
  sump = p1 + p2
  
  hellinger_distance = sqrt(1 - sqrt((2*sqrt(prodp))/(sump)) * exp(-(prodp*(m1 - m2)^2)/(4*(sump))))    
  return(as.numeric(hellinger_distance))
  
}

compute_hellinger_gammas <- function # Compute Gamma HD analitically
(hyperpar1, 
 ### Should be a named numeric array with 'shape' and 'rate' elements of 
 ### the first Gamma distribution.
 hyperpar2
 ### Should be a named numeric array with 'shape' and 'rate' elements of 
 ### the second Gamma distribution.
){
  
  a1 = hyperpar1[ ,"shape"]; b1 = hyperpar1[ ,"rate"]
  a2 = hyperpar2[ ,"shape"]; b2 = hyperpar2[ ,"rate"]
  
  hellinger_distance = sqrt(1 - gamma((a1 + a2)/2)*sqrt((b1^a1 * b2^a2)/(gamma(a1) * gamma(a2) * ((b1 + b2)/2)^(a1+a2))))
  return(as.numeric(hellinger_distance))
  
}

compute_hellinger_numerically <- function # Compute HD numerically
### Function used to compute the Hellinger distance when analytical
### form is not available or implemented.
(log_dist1,
 ### should be a function that takes 'x'  and 'hypepar' as arguments 
 ### and log_dist1(x, hyper) returns the log density value evaluated at 'x'. 
 hyperpar1,
 ### extra parameter to be passed to 'log_dist1', as hyperparameters for example.
 log_dist2,
 ### same as 'log_dist1'. 
 hyperpar2,
 ### same as 'hyperpar1'
 domain
 ### The end-points of the domain of 'log_dist1' and 'log_dist2', e.g. domain = c(-Inf, Inf). 
 ### In case the domain of 'log_dist1' and 'log_dist2' differ, include here the 'smaller' one, so that
 ### the integration is well defined.
 ){

  if (is.null(domain)) stop("'domain' is required for numerically computed HD.")

  g <- function(x, hyperpar1, hyperpar2){
    exp(0.5 *(log_dist1(x, hyperpar1) + log_dist2(x, hyperpar2))) 
  }
  bc_coef = integrate(f = g, lower = domain[1], upper = domain[2], 
                      rel.tol = getOption("priorSens_integrate_rel.tol"),
                      abs.tol = getOption("priorSens_integrate_abs.tol"),
                      subdivisions = getOption("priorSens_integrate_subdivisions"),
                      hyperpar1 = hyperpar1,
                      hyperpar2 = hyperpar2)$value
  bc_coef = min(bc_coef, 1)  
  
  hellinger_distance = sqrt(1 - bc_coef)
  return(as.numeric(hellinger_distance))
  
}

find_quantile <- function(quantile, log_dist, hyperpar, mean_dist = NULL){
  
  dist <- function(x, hyperpar){
    func = exp(log_dist(x, hyperpar))
    finite = is.finite(func)
    func[!finite] = 0
    return(func)
  }
  
  mean_func <- function(x, hyperpar) x*dist(x, hyperpar)
  
  CDF <- function(x) integrate(f=dist, lower=-Inf, upper=x, hyperpar = hyperpar)$value
  objective <- function(x, quantile) (CDF(x) - quantile)^2
  
  if (is.null(mean_dist)){
    
    mean_dist = integrate(f = mean_func, lower = -Inf, upper = Inf, hyperpar = hyperpar)$value  
    
  }
    
  result = nlminb(start=mean_dist, objective=objective, quantile = quantile)$par
  return(result)
  
}

compute_hellinger_numerically_grid_pertubation <- function # Compute HD numerically
### Function used to compute the Hellinger distance when analytical
### form is not available or implemented.
(log_dist,
 ### should be a function that takes 'x'  and 'hypepar' as arguments 
 ### and log_dist1(x, hyper) returns the log density value evaluated at 'x'. 
 hyperpar1,
 ### extra parameter to be passed to 'log_dist1', as hyperparameters for example.
 hyperpar2,
 ### same as 'hyperpar1'
 method = "integrate",
 ### Integration method used. Methods available are 'integrate' and 'simpson13rule'.
 domain = c(-Inf, Inf),
 ### Domain to be used witin the integration.
 number_points = 100
 ### Number of points using within "simpson13rule" integration method.
){
  
  
  if (method == "integrate"){
    
    g <- function(x, hyperpar1, hyperpar2){
      func = exp(0.5 *(log_dist(x, hyperpar1) + log_dist(x, hyperpar2)))
      finite = is.finite(func)
      func[!finite] = 0
      return(func)
    }
    
    bc_coef = integrate(f = g, lower = domain[1], upper = domain[2], 
                        rel.tol = getOption("priorSens_integrate_rel.tol"),
                        abs.tol = getOption("priorSens_integrate_abs.tol"),
                        subdivisions = getOption("priorSens_integrate_subdivisions"),
                        hyperpar1 = hyperpar1,
                        hyperpar2 = hyperpar2)$value
    
  } else if (method == "simpson13rule"){
    
    g <- function(x, hyperpar1, hyperpar2){
      exp(0.5 *(log_dist(x, hyperpar1) + log_dist(x, hyperpar2))) 
    }
    
    density_matrix = matrix(0, number_points, 2)
    density_matrix[,1] = seq(domain[1], domain[2], length.out=number_points)
    density_matrix[,2] = g(density_matrix[,1], hyperpar1, hyperpar2)
    
    bc_coef = simp13rule_matrix(density_matrix, apply_spline = FALSE)
    
  } else {
    
    stop("Wrong integration method.")
    
  }
  
  bc_coef = min(bc_coef, 1)  
  
  hellinger_distance = sqrt(1 - bc_coef)
  return(as.numeric(hellinger_distance))
  
} 
  
compute_hellinger_distance_inla <- function # Computation of the Hellinger distance between two inla posteriors
(inla_posterior_1, 
 ### inla posterior marginal 1
 inla_posterior_2,
 ### inla posterior marginal 2
 renormalize = TRUE,
 ### If TRUE (default), both 'inla_posterior_1' and 'inla_posterior_2'
 ### will renormalized to guarantee they integrate to 1.
 method = getOption("priorSens_inla_HD_integration_method"),
 ### integration method used to within the 'compute_hellinger_distance_inla'
 ### Currently, the following methods are available: 
 ### 'integrate', 'simpson13rule'. Default is 'simpson13rule'.
 integration_limits = TRUE
 ### If TRUE (default) use 10% rule for integration limits, 
 ### otherwise uses (-Inf, Inf).
){

  if (!all(inla_posterior_1[,1] == inla_posterior_2[,1])) 
    stop("compute_hellinger_distance_inla: Domains of the inla marginals should be the same.")
  
  if (renormalize){
    
    if (method == "integrate"){
    
    range_post_1 <- range(inla_posterior_1[,1])
    length_range_post_1 = range_post_1[2] - range_post_1[1]  
    
    cte1 = integrate(f=inla.dmarginal, 
                     range_post_1[1] - 0.1*length_range_post_1, 
                     range_post_1[2] + 0.1*length_range_post_1, 
                     rel.tol = getOption("priorSens_integrate_rel.tol"),
                     abs.tol = getOption("priorSens_integrate_abs.tol"),
                     subdivisions = getOption("priorSens_integrate_subdivisions"),
                     marginal = inla_posterior_1)$value

    range_post_2 <- range(inla_posterior_2[,1])
    length_range_post_2 = range_post_2[2] - range_post_2[1]  
    
    cte2 = integrate(f=inla.dmarginal, 
                     range_post_2[1] - 0.1*length_range_post_2, 
                     range_post_2[2] + 0.1*length_range_post_2,  
                     rel.tol = getOption("priorSens_integrate_rel.tol"),
                     abs.tol = getOption("priorSens_integrate_abs.tol"),
                     subdivisions = getOption("priorSens_integrate_subdivisions"),
                     marginal = inla_posterior_2)$value
    
    } else if (method == "simpson13rule"){
      
      cte1 = simp13rule_matrix(density_matrix = inla_posterior_1)
      cte2 = simp13rule_matrix(density_matrix = inla_posterior_2)
      
    } else {
      
      stop("Wrong integration method.")
      
    }
    
    inla_posterior_1[,2] = inla_posterior_1[,2]/cte1
    inla_posterior_2[,2] = inla_posterior_2[,2]/cte2
    ##details<<
    ## Before computing the hellinger distance, we renormalize both posteriors just
    ## to make sure they integrate to 1, this turns out to be inefficient if we
    ## already know they are properly normalized. In this case set 'renormalize = FALSE'
    
  }
  
  bc_integrand <- cbind(inla_posterior_1[,1], exp(0.5 * (log(inla_posterior_1[,2]) + log(inla_posterior_2[,2]))))
  
  if (method == "integrate"){
    
    range_bc_integrand <- range(bc_integrand[,1])
    length_range = range_bc_integrand[2] - range_bc_integrand[1]
    
    integration_limit_inf = - Inf
    integration_limit_sup = Inf
    
    if (integration_limits){
      integration_limit_inf = range_bc_integrand[1] - 0.1*length_range
      integration_limit_sup = range_bc_integrand[2] + 0.1*length_range
    }
    
    bc_coef <- integrate(f=inla.dmarginal, 
                         lower = integration_limit_inf, 
                         upper = integration_limit_sup, 
                         rel.tol = getOption("priorSens_integrate_rel.tol"),
                         abs.tol = getOption("priorSens_integrate_abs.tol"),
                         subdivisions = getOption("priorSens_integrate_subdivisions"),
                         marginal = bc_integrand)$value
    
  } else if (method == "simpson13rule"){
    
    bc_coef = simp13rule_matrix(density_matrix = bc_integrand)
    
  } else {
    
    stop("Wrong integration method.")
    
  }
  
  bc_coef = min(bc_coef, 1)  
  
  return(sqrt(1-bc_coef))
  ##value<<
  ## return the Hellinger distance.
  
}


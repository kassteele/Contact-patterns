\name{compute_roots_r_given_thetas}
\alias{compute_roots_r_given_thetas}
\title{Compute radial parameter values given a fixed theta.}
\description{This function returns radial parameters on the epsilon grid for a given theta parameter. }
\usage{compute_roots_r_given_thetas(new_theta, objective_function, hyperpar, 
    grid_epsilon, coordinate_correction = NULL, log_prior_density = NULL, 
    method = getOption("priorSens_root_method"), search_interval = c(-3, 
        3), integration_method = "integrate", domain = c(-Inf, 
        Inf), number_points = 100, ...)}
\arguments{
  \item{new_theta}{The theta parameter for which you want to find radial parameter values on the epsilon grid.
Needs to be between [-pi, pi].}
  \item{objective_function}{objective function used to compute the roots. Should take the form 
f(x, hyperpar, fixed_theta, grid_epsilon, coordinate_correction, log_prior_density) and
return (HD - grid_epsilon). There are also some pre-defined cases, like
objective_to_compute_roots_r_given_theta_gaussian for Gaussian, 
objective_to_compute_roots_r_given_theta_gamma for gamma and
objective_to_compute_roots_r_given_theta_custom for user-defined prior functions.}
  \item{hyperpar}{For the Gaussian case needs to be a named array with 'mean' and 'prec' elements.
For the Gamma case needs to be a named array with 'shape' and 'rate' elements.
For a user defined function should be consistent to what the defined function requires, 
but needs to be a numeric vector of length 2 defined on the real line.}
  \item{grid_epsilon}{epsilon used on the epsilon grid computations.}
  \item{coordinate_correction}{this will specify adjustments necessary to go from polar to cartesian coordinates
in order to cover the whole space efficiently. See function 'coord_correc_func' for details.}
  \item{log_prior_density}{log prior density function. For user defined prior distributions, it should be a 
function that takes 'x' defined on the 
real line and a numeric vector of length 2 called 'hypepar' with the hyperparameters of 
'log_prior_density' also defined on the real line as arguments.
log_prior_density(x, hypepar) must return the log density value evaluated at 'x'.
This argument should be NULL for particular cases already implemented, like 'gaussian'
and 'gamma'. }
  \item{method}{Methods used to find the radius of the polar coordinate system. 
Currently, the following methods are available: 'nlminb', 'uniroot.all', 'uniroot'. 
Default is 'uniroot'.}
  \item{search_interval}{Only if method equals 'uniroot.all' or 'uniroot'. A vector containing the end-points 
of the interval to be searched for the root of r in log scale. If no root in found, 
try to increase this interval. Remember that we have applied some kind of orthogonalization
in the polar coordinate system, so that r is expected to be around 1, so a 
search_interval = c(-3, 3) implies a search interval (0.049, 20.085) for r, which is likely
to contain the root.}
  \item{integration_method}{Integration method in case we use numerical integration.}
  \item{domain}{Necessary only when computing the HD numerically. The end-points of the hellinger 
distance integral.}
  \item{number_points}{Number of points used in the integration routine when simpson13rule is used.}
  \item{\dots}{Extra arguments for the 'method' used.}
}



\author{Thiago Guerrera Martins}






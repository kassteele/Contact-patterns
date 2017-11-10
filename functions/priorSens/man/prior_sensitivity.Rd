\name{prior_sensitivity}
\alias{prior_sensitivity}
\title{prior sensitivity}
\description{For Gaussian priors use 'gaussian', for gamma priors use 'gamma'.
For user defined prior distributions, it should be a 
function that takes 'x' defined on the 
real line and a numeric vector of length 2 called 'hypepar' with the hyperparameters of 
'log_prior_density' also defined on the real line as arguments.
log_prior_density(x, hypepar) must return the log density value evaluated at 'x'.}
\usage{prior_sensitivity(log_prior_density, hyperpar, inla_marginal_posterior, 
    number_axis_points = getOption("priorSens_number_axis_points"), 
    grid_epsilon = getOption("priorSens_grid_epsilon"), method = getOption("priorSens_root_method"), 
    integration_limits = TRUE, parallel = TRUE, ...)}
\arguments{
  \item{log_prior_density}{
}
  \item{hyperpar}{Depends on the objective function. For the Gaussian case needs to be a named 
array with 'mean' and 'prec' elements. For the Gamma case needs to be a named 
array with 'shape' and 'rate' elements.}
  \item{inla_marginal_posterior}{posterior marginal returned by INLA.}
  \item{number_axis_points}{Number of points used in the grid. Default is 400.}
  \item{grid_epsilon}{epsilon used in the definition of the contour. Default is 0.05.}
  \item{method}{Methods used to find the r. Currently, the following methods are available: 
'nlminb', 'uniroot.all', 'uniroot'. Default is 'uniroot'.}
  \item{integration_limits}{Set TRUE if you want to use 10\% rule for integration limits, otherwise uses Inf. Default is TRUE.}
  \item{parallel}{Use parallel computing through the multicore package.}
  \item{\dots}{Extra parameters for 'compute_roots_r_given_thetas' function.}
}



\author{Thiago Guerrera Martins}






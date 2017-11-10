\name{compute_grid_polar}
\alias{compute_grid_polar}
\title{Compute the epsilon grid using polar coordinates}

\usage{compute_grid_polar(number_axis_points, log_prior_density, hyperpar, 
    grid_epsilon = getOption("priorSens_grid_epsilon"), parallel = TRUE, 
    method = getOption("priorSens_root_method"), coord_explr_method = "nlminb", 
    coord_explr_method_options = NULL, integration_method = "integrate", 
    domain = c(-Inf, Inf), number_points = 100, ...)}
\arguments{
  \item{number_axis_points}{Number of points used in the grid.}
  \item{log_prior_density}{For Gaussian priors use 'gaussian', for gamma priors use 'gamma'.
For user defined prior distributions, it should be a 
function that takes 'x' defined on the 
real line and a numeric vector of length 2 called 'hypepar' with the hyperparameters of 
'log_prior_density' also defined on the real line as arguments.
log_prior_density(x, hypepar) must return the log density value evaluated at 'x'.}
  \item{hyperpar}{For the Gaussian case needs to be a named array with 'mean' and 'prec' elements.
For the Gamma case needs to be a named array with 'shape' and 'rate' elements.}
  \item{grid_epsilon}{epsilon used in the definition of the contour.}
  \item{parallel}{Use parallel computing through the multicore package.}
  \item{method}{Methods used to find the radius in the polar coordinate system. 
Currently, the following methods are available: 
'nlminb', 'uniroot.all', 'uniroot'. Default is 'uniroot'.}
  \item{coord_explr_method}{This is advanced option and sets the method used to find the root
under coordinate exploration mode, which is more delicate and 
therefore needs to use a more robust method.}
  \item{coord_explr_method_options}{A list with options to be passed to 'method_coord_exploration' routine.}
  \item{integration_method}{Integration method in case we use numerical integration.}
  \item{domain}{Necessary only when computing the HD numerically. The end-points of the hellinger 
distance integral.}
  \item{number_points}{Number of points used in the integration routine when simpson13rule is used.}
  \item{\dots}{Extra parameters for 'compute_roots_r_given_thetas' function.}
}



\author{Thiago Guerrera Martins}






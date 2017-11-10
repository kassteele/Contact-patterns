\name{objective_to_compute_roots_r_given_theta_custom}
\alias{objective_to_compute_roots_r_given_theta_custom}
\title{Return HD - epsilon between two user-defined priors}

\usage{objective_to_compute_roots_r_given_theta_custom(x, hyperpar, 
    fixed_theta, grid_epsilon, coordinate_correction, log_prior_density = NULL, 
    integration_method = "integrate", domain = c(-Inf, Inf), 
    number_points = 100)}
\arguments{
  \item{x}{x = log(r), where r is the radius.}
  \item{hyperpar}{should be a numeric vector, where the elements can take value on the real line.}
  \item{fixed_theta}{a given theta within [-pi, pi]}
  \item{grid_epsilon}{epsilon used on the epsilon grid computations.}
  \item{coordinate_correction}{this will specify adjustments necessary to go from polar to cartesian coordinates
in order to cover the whole space efficiently. See function 'coord_correc_func' for details.}
  \item{log_prior_density}{log prior density function. Should be a function that takes 'x' defined on the 
real line and a numeric vector of length 2 called 'hypepar' with the hyperparameters of 
'log_prior_density' also defined on the real line as arguments.
log_prior_density(x, hypepar) must return the log density value evaluated at 'x'.}
  \item{integration_method}{Integration method in case we use numerical integration.}
  \item{domain}{Necessary only when computing the HD numerically. The end-points of the hellinger 
distance integral.}
  \item{number_points}{Number of points used in the integration routine when simpson13rule is used.}
}



\author{Thiago Guerrera Martins}






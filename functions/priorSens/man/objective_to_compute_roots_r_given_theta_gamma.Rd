\name{objective_to_compute_roots_r_given_theta_gamma}
\alias{objective_to_compute_roots_r_given_theta_gamma}
\title{Return HD - epsilon between two Gamma priors}

\usage{objective_to_compute_roots_r_given_theta_gamma(x, hyperpar, fixed_theta, 
    grid_epsilon, coordinate_correction, log_prior_density = NULL, 
    ...)}
\arguments{
  \item{x}{x = log(r), where r is the radius.}
  \item{hyperpar}{should be a named array with 'shape' and 'rate' elements.}
  \item{fixed_theta}{a given theta within [-pi, pi]}
  \item{grid_epsilon}{epsilon used on the epsilon grid computations.}
  \item{coordinate_correction}{this will specify adjustments necessary to go from polar to cartesian coordinates
in order to cover the whole space efficiently. See function 'coord_correc_func' for details.}
  \item{log_prior_density}{log prior density function. Not necessary for the 'gamma' case.}
  \item{\dots}{
}
}



\author{Thiago Guerrera Martins}






\name{compute_sensitivity}
\alias{compute_sensitivity}
\title{Compute sensitivity measures for every point in the Grid}

\usage{compute_sensitivity(grid, inla_marginal_posterior, method = getOption("priorSens_correc_post_integration_method"), 
    method_hellinger = getOption("priorSens_hellinger_integration_method"), 
    integration_limits = TRUE)}
\arguments{
  \item{grid}{Grid values computed with 'compute_grid_polar' function.}
  \item{inla_marginal_posterior}{posterior marginal returned by INLA.}
  \item{method}{integration method used to compute the normalization cte of the corrected posteriors.}
  \item{method_hellinger}{integration method used to within the 'compute_hellinger_distance_inla'}
  \item{integration_limits}{Set TRUE if you want to use 10\% rule for integration limits, otherwise uses Inf. Default is TRUE.}
}



\author{Thiago Guerrera Martins}






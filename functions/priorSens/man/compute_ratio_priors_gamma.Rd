\name{compute_ratio_priors_gamma}
\alias{compute_ratio_priors_gamma}
\title{Compute the ratio 'new_prior/base_prior' for gamma case.}

\usage{compute_ratio_priors_gamma(hyperpar_new, hyperpar_old, inla_internal_marginal)}
\arguments{
  \item{hyperpar_new}{Hyperparameter of the new prior. Needs to be a named numerical array with 'shape' 
and 'rate' parameters.}
  \item{hyperpar_old}{Hyperparameter of the old/base prior. Needs to be a named numerical array with 'shape' 
and 'rate' parameters.}
  \item{inla_internal_marginal}{Needs to be the inla posterior marginal on the theta scale where theta is assigned a 
Gamma prior.}
}



\author{Thiago Guerrera Martins}






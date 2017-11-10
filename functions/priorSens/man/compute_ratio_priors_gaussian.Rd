\name{compute_ratio_priors_gaussian}
\alias{compute_ratio_priors_gaussian}
\title{Compute the ratio 'new_prior/base_prior' for Gaussian case.}

\usage{compute_ratio_priors_gaussian(hyperpar_new, hyperpar_old, inla_internal_marginal)}
\arguments{
  \item{hyperpar_new}{Hyperparameter of the new prior. Needs to be a named numerical array with 'mean' and 
'prec' parameters.}
  \item{hyperpar_old}{Hyperparameter of the old/base prior. Needs to be a named numerical array with 'mean' and 
'prec' parameters.}
  \item{inla_internal_marginal}{Needs to be the inla posterior marginal on the theta scale where theta is assigned a 
Gaussian prior.}
}



\author{Thiago Guerrera Martins}






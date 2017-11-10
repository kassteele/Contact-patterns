\name{compute_corrected_posteriors}
\alias{compute_corrected_posteriors}
\title{Compute corrected posteriors from INLA marginals.}
\description{This function computes corrected posteriors for hyperparameters values of the epsilon-grid given
an INLA posterior marginal in the appropriate representation.}
\usage{compute_corrected_posteriors(prior_type, hyperpar_new, hyperpar_old, 
    inla_marginal_posterior, method = getOption("priorSens_correc_post_integration_method"), 
    integration_limits = TRUE)}
\arguments{
  \item{prior_type}{A character vector indication what is the family of the prior. Currently, available for 
"gaussian" and "gamma".  }
  \item{hyperpar_new}{hyperparameters of the new prior, see Details. }
  \item{hyperpar_old}{hyperparameters of the old/base prior, see Details. }
  \item{inla_marginal_posterior}{Needs to be the inla posterior marginal on the theta scale where theta is assigned a 
'prior_type'.}
  \item{method}{Integration method used to normalize the corrected posterior}
  \item{integration_limits}{Set TRUE if you want to use 10\% rule for integration limits, otherwise uses Inf. Default is TRUE.}
}
\details{Currently, this function is implemented only for 'prior_type' equal to 'gaussian' and 'gamma'. 

If prior_type="gaussian", this includes a named numerical array with 'mean' and 'prec' 
parameters. If prior_type="gamma", this includes a named numerical array with
'shape' and 'rate' parameters.}
\value{Corrected posterior marginals using INLA format.}

\author{Thiago Guerrera Martins}






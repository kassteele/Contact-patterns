\name{compute_hellinger_distance_inla}
\alias{compute_hellinger_distance_inla}
\title{Computation of the Hellinger distance between two inla posteriors}

\usage{compute_hellinger_distance_inla(inla_posterior_1, inla_posterior_2, 
    renormalize = TRUE, method = "integrate", integration_limits = TRUE)}
\arguments{
  \item{inla_posterior_1}{inla posterior marginal 1}
  \item{inla_posterior_2}{inla posterior marginal 2}
  \item{renormalize}{If TRUE (default), both 'inla_posterior_1' and 'inla_posterior_2'
will renormalized to guarantee they integrate to 1.}
  \item{method}{integration method used.}
  \item{integration_limits}{If TRUE (default) use 10\% rule for integration limits, 
otherwise uses (-Inf, Inf).}
}
\details{Before computing the hellinger distance, we renormalize both posteriors just
to make sure they integrate to 1, this turns out to be inefficient if we
already know they are properly normalized. In this case set 'renormalize = FALSE'}
\value{return the Hellinger distance.}

\author{Thiago Guerrera Martins}






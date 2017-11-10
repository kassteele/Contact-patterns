\name{compute_hellinger_numerically}
\alias{compute_hellinger_numerically}
\title{Compute HD numerically}
\description{Function used to compute the Hellinger distance when analytical
form is not available or implemented.}
\usage{compute_hellinger_numerically(log_dist1, hyperpar1, log_dist2, 
    hyperpar2, domain)}
\arguments{
  \item{log_dist1}{should be a function that takes 'x'  and 'hypepar' as arguments 
and log_dist1(x, hyper) returns the log density value evaluated at 'x'. }
  \item{hyperpar1}{extra parameter to be passed to 'log_dist1', as hyperparameters for example.}
  \item{log_dist2}{same as 'log_dist1'. }
  \item{hyperpar2}{same as 'hyperpar1'}
  \item{domain}{The end-points of the domain of 'log_dist1' and 'log_dist2', e.g. domain = c(-Inf, Inf). 
In case the domain of 'log_dist1' and 'log_dist2' differ, include here the 'smaller' one, so that
the integration is well defined.}
}



\author{Thiago Guerrera Martins}






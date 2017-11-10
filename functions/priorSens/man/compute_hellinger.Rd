\name{compute_hellinger}
\alias{compute_hellinger}
\title{Hellinger distance between two distributions}
\description{Computes the hellinger distance between two distributions numerically.
The computation is made analitically in case both distributions are Gaussians
or Gammas.}
\usage{compute_hellinger(log_dist1, hyperpar1, log_dist2, hyperpar2, 
    integration_method = "integrate", domain = c(-Inf, Inf), 
    number_points = 100)}
\arguments{
  \item{log_dist1}{should be a function that takes 'x' defined in the real line and 'hypepar' as arguments 
and log_dist1(x, hyper) returns the log density value evaluated at 'x'. Write "gaussian" 
or "gamma" in case both 'dist1' and 'dist2' are Gaussian or gamma distributions respectively.}
  \item{hyperpar1}{extra parameter to be passed to 'log_dist1', as hyperparameters for example. 
If "gaussian", this includes a named numerical array with 'mean' and 'prec' 
parameters. If "gamma", this includes a named numerical array with
'shape' and 'rate' parameters.}
  \item{log_dist2}{same as 'log_dist1'. }
  \item{hyperpar2}{same as 'hyperpar1'}
  \item{integration_method}{Integration method in case we use numerical integration.}
  \item{domain}{Necessary only when computing the HD numerically. The end-points of the hellinger 
distance integral.}
  \item{number_points}{Number of points used in the integration routine when simpson13rule is used.}
}



\author{Thiago Guerrera Martins}






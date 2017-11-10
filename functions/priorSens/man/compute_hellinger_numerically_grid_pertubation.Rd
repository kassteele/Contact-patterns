\name{compute_hellinger_numerically_grid_pertubation}
\alias{compute_hellinger_numerically_grid_pertubation}
\title{Compute HD numerically}
\description{Function used to compute the Hellinger distance when analytical
form is not available or implemented.}
\usage{compute_hellinger_numerically_grid_pertubation(log_dist, hyperpar1, 
    hyperpar2, method = "integrate", domain = c(-Inf, Inf), number_points = 100)}
\arguments{
  \item{log_dist}{should be a function that takes 'x'  and 'hypepar' as arguments 
and log_dist1(x, hyper) returns the log density value evaluated at 'x'. }
  \item{hyperpar1}{extra parameter to be passed to 'log_dist1', as hyperparameters for example.}
  \item{hyperpar2}{same as 'hyperpar1'}
  \item{method}{Integration method used.}
  \item{domain}{Domain to be used witin the integration.}
  \item{number_points}{Number of points using within "simpson13rule" integration method.}
}



\author{Thiago Guerrera Martins}






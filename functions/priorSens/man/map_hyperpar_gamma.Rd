\name{map_hyperpar_gamma}
\alias{map_hyperpar_gamma}
\title{Map from polar coordinate to cartesian coordinate for Gamma prior.}

\usage{map_hyperpar_gamma(hyperpar, r, theta, coord_correc_x, coord_correc_y)}
\arguments{
  \item{hyperpar}{The initial hyperparameters of the Gamma prior. Needs to be a numeric 
vector where the first element is the shape and the second element is the rate.}
  \item{r}{radius of the polar coordinate.}
  \item{theta}{The degree of the polar coordinate.}
  \item{coord_correc_x}{x-axis correction, see function 'coord_correc_func' for details.}
  \item{coord_correc_y}{y-axis correction, see function 'coord_correc_func' for details.}
}



\author{Thiago Guerrera Martins}






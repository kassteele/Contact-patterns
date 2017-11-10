\name{map_hyperpar_custom}
\alias{map_hyperpar_custom}
\title{Map from polar coordinate to cartesian coordinate for an user-defined prior.}

\usage{map_hyperpar_custom(hyperpar, r, theta, coord_correc_x, coord_correc_y)}
\arguments{
  \item{hyperpar}{The initial hyperparameters, defined on the real line, of the user-defined prior. Needs to be a numeric 
vector of length 2.}
  \item{r}{radius of the polar coordinate.}
  \item{theta}{The degree of the polar coordinate.}
  \item{coord_correc_x}{x-axis correction, see function 'coord_correc_func' for details.}
  \item{coord_correc_y}{y-axis correction, see function 'coord_correc_func' for details.}
}



\author{Thiago Guerrera Martins}






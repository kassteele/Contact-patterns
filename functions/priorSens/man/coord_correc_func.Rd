\name{coord_correc_func}
\alias{coord_correc_func}
\title{Specify which coordinate correction to use.}

\usage{coord_correc_func(theta, coordinate_correction)}
\arguments{
  \item{theta}{theta parameter of the polar coordinate, within [-pi, pi].}
  \item{coordinate_correction}{NULL if no coordinate correction is to be applied. Otherwise, a list
with the following elements: 'xneg', 'xpos', 'yneg', 'ypos' indicating
the coordinate corrections that will be applied dependent on the value of theta, 
see Details.}
}
\details{The theta space is divided into 4 quadrants: if theta belongs to 
[-pi, -pi/2] then we are in the negative side of the x-xis and of the y-axis.
if theta belongs to [-pi/2, 0] we are on the positive side of the x-axis and on
the negative side of the y axis, and so on.}


\author{Thiago Guerrera Martins}






% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geometryfns.R
\name{real_line}
\alias{real_line}
\title{real line}
\usage{
real_line(measure = Euclid_dist(dim = 1L))
}
\arguments{
\item{measure}{an object of class \code{measure}}
}
\description{
Initialisation of the real-line (1D) manifold.
}
\details{
A real line is initialised using a \code{measure} object. By default, the measure object (\code{measure}) describes the distance between two points as the absolute difference between the two coordinates.
}
\examples{
R <- real_line()
print(type(R))
print(sp::dimensions(R))
}

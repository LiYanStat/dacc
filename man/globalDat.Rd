% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{globalDat}
\alias{globalDat}
\title{Sample Dataset Used in Numerical Studies of 
"Detection and Attribution Analysis of Temperature Changes with Estimating Equations".}
\format{
A data list with the observed and simulated data on global scale.
}
\usage{
data(globalDat)
}
\description{
A list of the observations and expected responses 
to different external forcings with name
\code{Y}, \code{X}, \code{ctlruns}, \code{nruns.X} and \code{Xtilde}
where
\itemize{
    \item \code{Y}: the \eqn{5^\circ \times 5^\circ} gridded observations on global scale
    \item \code{X}: a data matrix of the expected responses to external forcing;
    \item \code{ctlruns}: replicates of control runs from pre-industrial simulations
    \item \code{nruns.X}: number of runs for the estimated responses to external forcings
    \item \code{Xtilde}: the selected estimated responses to external forcing ANT and NAT
}
}
\examples{
data(globalDat)
}
\keyword{datasets}

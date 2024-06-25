##' Sample Dataset in Simulation Studies of 
##' "Regularized Fingerprinting in Detection and Attribution of Climate Change 
##' with Weight Matrix Optimizing the Efficiency in Scaling Factor Estimation".
##'
##' A data list of designed covariance matrix and the expected responses 
##' to the two forcings ANT and NAT with name
##' \code{Cov} and \code{X}
##' where
##' \itemize{
##'     \item \code{Cov}: a list of the true covariance matrices;
##'     \item \code{X}: a data matrix of the expected responses to external forcing;
##' }
##'
##' @docType data
##' @name simDat
##' @usage data(simDat)
##' @format A data list with two separate data sets.
##' @keywords datasets
##' @examples
##' data(simDat)
"simDat"


##' Sample Dataset Used in Numerical Studies of 
##' "Detection and Attribution Analysis of Temperature Changes with Estimating Equations".
##'
##' A list of the observations and expected responses 
##' to different external forcings with name
##' \code{Y}, \code{X}, \code{ctlruns}, \code{nruns.X} and \code{Xtilde}
##' where
##' \itemize{
##'     \item \code{Y}: the \eqn{5^\circ \times 5^\circ} gridded observations on global scale
##'     \item \code{X}: a data matrix of the expected responses to external forcing;
##'     \item \code{ctlruns}: replicates of control runs from pre-industrial simulations
##'     \item \code{nruns.X}: number of runs for the estimated responses to external forcings
##'     \item \code{Xtilde}: the selected estimated responses to external forcing ANT and NAT
##' }
##'
##' @docType data
##' @name globalDat
##' @usage data(globalDat)
##' @format A data list with the observed and simulated data on global scale.
##' @keywords datasets
##' @examples
##' data(globalDat)
"globalDat"
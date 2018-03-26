##' 5-year mean Rx1day precipitation over four large regions in Canada. 
##' Including obvered data from HadEx2 and simulated data from CMIP5 model.
##'
##' A data list of simulated and observed data with elements named
##' \code{x.all}, \code{x.ant}, \code{x.nat} and \code{Y.observation}, 
##' where
##' \itemize{
##'     \item \code{x.all}: simulated replicates under ALL forcing;
##'     \item \code{x.ant}: simulated replicates under ANT forcing;
##'     \item \code{x.nat}: simulated replicates under NAT forcing;
##'     \item \code{Y.observation}: observed precipitation data.
##' }
##' 
##' @docType data
##' @name 
##' @usage data(CanadaPrcp)
##' @format A data list with four separate data sets.
##' @keywords datasets
##' @examples
##' data(CanadaPrcp)
"CanadaPrcp"

##' Sample Data of Observed and Simulated Global Mean Temperature.
##'
##' A data list of simulated and observed data with elements named
##' \code{ctlruns} and \code{obs}, \code{signal},
##' where
##' \itemize{
##'     \item \code{ctlruns}: a matrix of the simulated 
##'     global mean temperature;
##'     \item \code{obs}: a vector of the observed data;
##'     \item \code{signal}: a matrix of the signal factors.
##' }
##'
##' @docType data
##' @name cnrmDat
##' @usage data(cnrmDat)
##' @format A data list with three separate data sets.
##' @keywords datasets
##' @examples
##' data(cnrmDat)
"cnrmDat"

##' Sample Data of Observed and Simulated Global Mean Temperature under T4
##' truncation.
##'
##' A data list of simulated and observed data with elements named
##' \code{ctlruns} and \code{obs}, \code{signal},
##' where
##' \itemize{
##'     \item \code{ctlruns}: a matrix of the simulated 
##'     global mean temperature;
##'     \item \code{obs}: a vector of the observed data;
##'     \item \code{signal}: a matrix of the signal factors.
##' }
##' 
##' @docType data
##' @name HistT4Dat
##' @usage data(HistT4Dat)
##' @format A data list with three separate data sets.
##' @keywords datasets
##' @examples
##' data(HistT4Dat)
"HistT4Dat"
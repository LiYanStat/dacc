##' Optimal Fingerprinting via total least square regression.
##'
##' This function estimates the signal factors and corresponding confidence 
##' interval via the estimating equation or total least squares.
##'
##' @param Xtilde \eqn{n \times p} matrix, signal pattern to be detected.
##' @param Y \eqn{n \times 1} matrix, length \eqn{S \times T}, observed climate variable.
##' @param mruns number of ensembles to estimate the corresponding pattern. 
##' It is used as the scale of the covariance matrix for \eqn{X_i}.
##' @param ctlruns.sigma \eqn{m \times n} matrix, a group of \eqn{m} independent control 
##' runs for estimating covariance matrix, which is used in point estimation of 
##' the signal factors.
##' @param ctlruns.bhvar \eqn{m \times n} matrix, another group of \eqn{m} independent control 
##' runs for estimating the corresponding confidence interval of the signal factors, 
##' in EE or PBC approach should be same as ctlruns.sigma.
##' @param S number of locations for the observed responses.
##' @param T number of time steps for the observed responses.
##' @param B number of replicates in bootstrap procedure, mainly for the PBC and TS methods, can be 
##' specified in "EE" method but not necessary. By default B = 0 as the default method is "EE".
##' @param Proj The projection matrix for computing for scaling factors of other external forcings 
##' with the current input when using EE. For example, when ALL and NAT are used for modeling, 
##' specifying the Proj matrix to return the results for ANT and NAT.
##' @param method for estimating the scaling factors and corresponding confidence interval
##' @param cov.method method for estimation of covariance matrix in confidence 
##' interval estimation of PBC method. (only for PBC method).
##' @param conf.level confidence level for confidence interval estimation.
##' @param missing indicator for whether missing values present in Y.
##' @param cal.a indicator for calculating the a value, otherwise use default value a = 1. (only for EE method)
##' @param ridge shrinkage value for adjusting the method for missing observations if missing = TRUE. (only for EE method)
##' @return a list of the fitted model including point estimate and
##' interval estimate of coefficients and corresponding estimate of standard error.
##' @author Yan Li
##' @keywords regression tls fingerprinting estimating equation
##' @references \itemize{ 
##' \item  Gleser (1981), Estimation in a Multivariate "Errors in Variables" 
##' Regression Model: Large Sample Results, \emph{Ann. Stat.} 9(1) 24--44.
##' \item Golub and Laon (1980), An Analysis of the Total Least Squares Problem,
##' \emph{SIAM J. Numer. Anal}. 17(6) 883--893.
##' \item Pesta (2012), Total least squares and bootstrapping with 
##' applications in calibration, \emph{Statistics} 47(5), 966--991.
##' \item Li et al (2021), Uncertainty in Optimal Fingerprinting is Underestimated, 
##'       \emph{Environ. Res. Lett.} 16(8) 084043.
##' \item Sai et al (2023), Optimal Fingerprinting with Estimating Equations,
##'       \emph{Journal of Climate} 36(20), 7109–-7122.
##' \item Li et al (2024), Detection and Attribution Analysis of Temperature Changes with Estimating Equations,
##'       \emph{Submitted to Journal of Climate}.
##' }
##' @examples
##' ## load the example dataset
##' data(simDat)
##' Cov <- simDat$Cov[[1]]
##' ANT <- simDat$X[, 1]
##' NAT <- simDat$X[, 2]
##' 
##' ## generate the simulated data set
##' ## generate regression observation
##' Y <- MASS::mvrnorm(n = 1, mu = ANT + NAT, Sigma = Cov)
##' ## generate the forcing responses
##' mruns <- c(1, 1)
##' Xtilde <- cbind(MASS::mvrnorm(n = 1, mu = ANT, Sigma = Cov / mruns[1]),
##'                MASS::mvrnorm(n = 1, mu = NAT, Sigma = Cov / mruns[2]))
##' ## control runs
##' ctlruns <- MASS::mvrnorm(100, mu = rep(0, nrow(Cov)), Sigma = Cov)
##' ## ctlruns.sigma for the point estimation and ctlruns.bhvar for the interval estimation
##' ctlruns.sigma <- ctlruns.bhvar <- ctlruns
##' ## number of locations
##' S <- 25
##' ## number of year steps
##' T <- 10
##' 
##' ## call the function to estimate the signal factors via EE
##' fingerprint(Xtilde, Y, mruns,
##'           ctlruns.sigma, ctlruns.bhvar,
##'           S, T,
##'           ## B = 0, by default
##'           method = "EE",
##'           conf.level = 0.9,
##'           cal.a = TRUE,
##'           missing = FALSE, ridge = 0)
##' @import magrittr
##' @importFrom MASS mvrnorm
##' @importFrom stats cov qnorm quantile sd var pnorm
##' @importFrom utils tail
##' @importFrom pracma ceil
##' @importFrom janitor remove_empty
##' @importFrom grDevices boxplot.stats
##' @export fingerprint

fingerprint <- function(Xtilde, Y, mruns, 
                        ctlruns.sigma, ctlruns.bhvar,
                        S, T, B = 0,
                        Proj = diag(ncol(Xtilde)),
                        method = c("EE", "PBC", "TS"),
                        cov.method = c("l2", "mv"),
                        conf.level = 0.90,
                        missing = FALSE,
                        cal.a = TRUE, ridge = 0) {
  Call <- match.call()
  Xtilde <- as.matrix(Xtilde)
  method <- match.arg(method)  ## method for the estimating procedure
  cov.method <- match.arg(cov.method)  ## method for the covariance matrix extimation in PBC and TS
  if(! method %in% c("EE", "PBC", "TS")) {
    stop("Unknow method for signal estimation")
  }
  if(is.null(colnames(Xtilde))) {
    colnames(Xtilde) <- paste0("forcings ", 1:ncol(Xtilde))
  }
  ## for the EE method
  if(method == "EE") {
    ## check the missing value
    if(any(is.na(Y))) {
      if(!missing) {
        message("missing values in observations")
      }
      missing <- TRUE
    } else {
      if(missing) {
        message("no missing values in observations")
      }
      missing <- FALSE
    }
    output <- fingerprintCEE(
      Xtilde = Xtilde, Y = Y, mruns = mruns,
      ctlruns.1 = ctlruns.sigma, 
      ctlruns.2 = ctlruns.bhvar,
      Proj = Proj,
      nS = S, nT = T, nB = B,
      conf.level = conf.level,
      cal.a = cal.a, missing = missing,
      ridge = ridge
    )
  } else {
    if (cov.method == "l2") {
      cov <- Covest(ctlruns.sigma, method = "l2")$output
    } else if (cov.method == "mv") {
      cov <- Covest(ctlruns.sigma, method = "mv", bandwidth = 0.35)$output
    } else {
      stop("Unknow method for covariance matrix estimation")
    }
    if(missing) {
      ## check the missing value
      nomis <- which(!is.na(Y))
    }
    ## for the PBC method
    if(method == "PBC") {
      output <- fingerprintTLS(
        X = Xtilde[nomis, , drop = FALSE], 
        Y = Y[nomis, , drop = FALSE],
        cov = cov[nomis, nomis, drop = FALSE], nruns.X = mruns,
        ctlruns = ctlruns.sigma[, nomis, drop = FALSE],
        precision = FALSE,
        conf.level = conf.level,
        conf.method = "PBC",  ## only output results from PBC (TSB also available)
        cov.method = cov.method,
        B = B)
    }
    ## for the two sample approach
    if(method == "TS") {
      output <- fingerprintTS(X = Xtilde[nomis, , drop = FALSE],
                              Y = Y[nomis, , drop = FALSE],
                              nruns.X = mruns, cov = cov[nomis, nomis, drop = FALSE],
                              Z.2 = ctlruns.bhvar[, nomis, drop = FALSE],
                              precision = FALSE,
                              conf.level = conf.level,
                              B = B)
    }
  }
  ## return the output
  return(output)
}
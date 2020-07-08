##' Optimal Fingerprinting via total least square regression.
##'
##' This function estimates the ratio of underestimation when using 
##' regularized covariance matrix as weight matrix in GTLS regression. 
##' The estimation of scaling ratio is computed from conduct the same 
##' estimation procedure to simulated replicates from the estimated 
##' regression model and get the ratio between empirical standard errors 
##' and formula based estimation of standard errors. 
##'
##' @param X signal pattern to be detected.
##' @param Y observed data.
##' @param nruns.X number of ensembles to estimate the corresponding pattern. 
##' It is used as the scale of the covariance matrix for Xi.
##' @param cov Weight matrix used in prewhitening process, which should be 
##' estimate of covariance matrix.
##' @param rep.num the sample size of independent control runs for estimating 
##' covariance matrix.
##' @param method method for estimation of covariance matrix in calibration 
##' method. It should be consistent to the method to get cov. 
##' @param B number of replicates in calibration, default value 
##' is 1000.
##' @return a list of the fitted model including point estimate and
##' interval estimate of coefficients and corresponding estimate of 
##' standard error.
##' @author Yan Li
##' @references \itemize{ 
##' \item  Gleser, Estimation in a Multivariate "Errors in Variables" 
##' Regression Model: Large Sample Results, 1981, Ann. Stat.
##' \item Golub and Laon, An Analysis of the Total Least Squares Problem,
##' 1980, SIAM J. Numer. Anal.
##' \item Pesta, Total least squares and bootstrapping with 
##' applications in calibration, 2012, Statistics.}
##' @examples
##' data(simDat)
##' ## set the true covariance matrix and expected pattern
##' Cov <- simDat$Cov[[1]]
##' ANT <- simDat$X[, 1]
##' NAT <- simDat$X[, 2]
##' ## estimate the covariance matrix
##' Z <- MASS::mvrnorm(100, mu = rep(0, nrow(Cov)), Sigma = Cov)
##' ## linear shrinkage estimator under l2 loss
##' Cov.est <- Covest(Z, method = "l2")$output
##' ## generate regression observation and pattern
##' nruns.X <- c(1, 1)
##' Y <- MASS::mvrnorm(n = 1, mu = ANT + NAT, Sigma = Cov)
##' X <- cbind(MASS::mvrnorm(n = 1, mu = ANT, Sigma = Cov),
##'            MASS::mvrnorm(n = 1, mu = NAT, Sigma = Cov))
##' res <- calibFingerprint(X, Y, nruns.X, rep.num = nrow(Z), Cov.est, method = "l2", B = 10)
##' @import MASS stats utils methods
##' @export calibFingerprint

calibFingerprint <- function(X, Y, nruns.X, rep.num, cov, method, B = 1000) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  ## check the colnames
  if(is.null(colnames(X))) {
    colnames(X) <- paste0("forcings ", 1:ncol(X))
  }
  if(is.null(colnames(Y))) {
    colnames(Y) <- "response"
  }
  tmpMat <- eigen(cov)
  cov.sinv <- Re(tmpMat$vectors %*% diag(1 / sqrt(tmpMat$values)) %*% 
                   t(tmpMat$vectors))
  
  tls <- function(X, Y, nruns.X) {
    if (! is.matrix(Y)) {
      stop("Y should be a n*1 matrix")
    }
    if (dim(X)[1] != dim(Y)[1])  {  ## check size of X and Y
      stop("sizes of inputs X, Y are not consistent")
    }
    n <- dim(Y)[1]  ## number of observations
    m <- dim(X)[2]  ## number of predictors
    ## Normalise the variance of X
    X <- X * t(sqrt(nruns.X) %*% matrix(1, 1, n))
    Estls <- function(X, Y) {
      M <- cbind(X, Y)
      svd.M <- svd(M)
      output <- svd.M$u %*% diag(c(svd.M$d[1:m], 0)) %*% t(svd.M$v)
      colnames(output) <- c(colnames(X), colnames(Y))
      output
    }
    output <- Estls(X, Y)
    output[, colnames(X)] <- output[, colnames(X)] / 
      t(sqrt(nruns.X) %*% matrix(1, 1, n))
    output
  }
  
  output <- tls(cov.sinv %*% X, cov.sinv %*% Y, nruns.X)
  output <- Re(tmpMat$vectors %*% diag(sqrt(tmpMat$values)) %*% 
                 t(tmpMat$vectors)) %*% output
  
  out.beta <- out.var <- NULL
  for (i in 1:B) {
    ## generate new dataset
    for(er in 1:500) {
      tmp.new <- tryCatch({
        Y.new <- mvrnorm(n = 1, mu = output[, colnames(Y)], Sigma = cov)
        X.new <- NULL
        for(X.ind in 1:ncol(X)) {
          X.new <- cbind(X.new,
                         mvrnorm(n = 1, mu = output[, X.ind], Sigma = cov / nruns.X[X.ind]))
        }
        ctlruns.new <- genClt(rep.num, 1, cov)[[1]]
        if (method == "l2") {
          Cov.new <- Covest(ctlruns.new, method = "l2")$output
        } else if (method == "mv") {
          Cov.new <- Covest(ctlruns.new, method = "mv", bandwidth = 0.35)$output
        }
        tlsFingerprint(X.new, Y.new, Cov.new, nruns.X, ctlruns.new, method = method, tsb = FALSE)
      }, error = function(e) {
        ""
      })
      if (tmp.new[1] != "") {
        break
      }
    }
    out.beta <- rbind(out.beta, tmp.new$coefficient)
    out.var <- rbind(out.var, as.vector(tmp.new$var.est))
  }
  list(rbind(apply(out.beta, 2, sd), apply(out.beta, 2, sd)) / matrix(colMeans(out.var, na.rm = TRUE), nrow = 2), 
       out.beta, out.var)
}

## generate independent replicates
genClt <- function(n, B, Cov) {
  lapply(1:B, 
         function(x) {
           MASS::mvrnorm(n, mu = rep(0, nrow(Cov)), Sigma = Cov)
         })
}

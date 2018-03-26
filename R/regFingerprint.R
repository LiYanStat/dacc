##' Optimal Fingerprinting via linear regression.
##'
##' This function detects the signal factors on the observed data via linear
##' regression.
##' 
##' @param X signal pattern to be detected.
##' @param Y observed data.
##' @param weights Weight matrix used in GLS, can be estimate of 
##' covariance matrix or precision matrix.
##' @param conf.level confidence level for computing confidence interval.
##' @param precision indicator for precision matrix, if precision 
##' matrix estimate is used, precision should be set to TRUE.
##' @return a list of the fitted model including point estimate and
##' interval estimate of coefficients and corresponding estimate of 
##' covariance matrix.
##' @author Yan Li
##' @keywords Regularized Fingerprinting, GLS
##' @examples
##' data(CanadaPrcp)
##' Y <- CanadaPrcp$Y.observation
##' X <- colMeans(CanadaPrcp$x.ant)
##' regFingerprint(X, Y, weights = diag(length(Y)))
##' @import MASS stats utils methods
##' @importFrom expm sqrtm
##' @export regFingerprint
regFingerprint <- function(X, Y, weights, conf.level = 0.95, precision = FALSE) {
  X <- as.matrix(X)
  if(! precision) {
    Df <- length(Y) - ncol(X)
    beta.hat <- solve(t(X) %*% solve(weights) %*% X) %*% t(X) %*% solve(weights) %*% Y
    ## covariance matrix
    ## cov.betahat <- solve(t(X) %*% solve(weights) %*% X) %*% 
    ##    t(X) %*% solve(weights) %*% Cov %*% solve(weights) %*% X %*% 
    ##    solve(t(X) %*% solve(weights) %*% X)
    sig.hat  <- as.numeric((t(Y - X %*% beta.hat) %*% solve(weights) %*% 
                              (Y - X %*% beta.hat) / Df))
    cov.betahat <-  sig.hat * (solve(t(X) %*% solve(weights) %*% X))
    CI.hat <- ConfInt(beta.hat, cov.betahat, conf.level = conf.level, Df = Df)
    
    list(beta.hat = beta.hat, cov.betahat = cov.betahat, CI.hat = CI.hat)
  } else {
    Df <- length(Y) - ncol(X)
    beta.hat <- solve(t(X) %*% weights %*% X) %*% t(X) %*% weights %*% Y
    sig.hat  <- as.numeric((t(Y - X %*% beta.hat) %*% weights %*% 
                              (Y - X %*% beta.hat) / Df))
    cov.betahat <-  sig.hat * (solve(t(X) %*% weights %*% X))
    CI.hat <- ConfInt(beta.hat, cov.betahat, conf.level = conf.level, Df = Df)
    list(beta.hat = beta.hat, cov.betahat = cov.betahat, CI.hat = CI.hat)
  }
}

## compute the confidence interval via student t distribution
ConfInt <- function(param, cov, conf.level, Df) {
  var <- sqrt(diag(cov))
  upper <- param + qt(p = conf.level + (1 - conf.level) / 2, df = Df) * var
  lower <- param - qt(p = conf.level + (1 - conf.level) / 2, df = Df) * var
  cbind(lower, upper)
}

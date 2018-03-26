##' Optimal Fingerprinting via total least square regression.
##'
##' This function detects the signal factors on the observed data via total 
##' least square linear regression model.
##'
##' @param X signal pattern to be detected.
##' @param Y observed data.
##' @param cov Weight matrix used in prewhitening process, can be estimate of 
##' covariance matrix or precision matrix.
##' @param nruns.X number of ensembles to estimate the corresponding pattern. 
##' It is used as the scale of the covariance matrix for Xi.
##' @param precision indicator for precision matrix, if precision 
##' matrix estimate is used, precision should be set to TRUE.
##' @return a list of the fitted model including point estimate and
##' interval estimate of coefficients and corresponding estimate of 
##' standard error.
##' @author Yan Li
##' @keywords Regularized Fingerprinting, Generalized TLS
##' @references \itemize{ 
##' \item  Gleser, Estimation in a Multivariate "Errors in Variables" 
##' Regression Model: Large Sample Results, 1981, Ann. Stat.
##' \item Golub and Laon, An Analysis of the Total Least Squares Problem,
##' 1980, SIAM J. Numer. Anal.
##' \item Pesta, Total least squares and bootstrapping with 
##' applications in calibration, 2012, Statistics.}
##' @examples
##' data(CanadaPrcp)
##' ## estimate covariance matrix
##' Z <- scale(CanadaPrcp$x.all, center = TRUE, scale = FALSE)
##' ## linear shrinkage estimator under l2 loss
##' Cov.est <- Covest(Z, method = "l2")$output
##' ## nonlinear shrinkage estimator under minimum variance loss
##' Cov.est <- Covest(Z, method = "mv")$output
##' ## regression obervation and pattern
##' Y <- CanadaPrcp$Y.observation
##' X <- colMeans(CanadaPrcp$x.ant)
##' tlsFingerprint(X, Y, nruns.X = length(X), Cov.est)
##' @import MASS stats utils methods
##' @importFrom expm sqrtm
##' @export tlsFingerprint
tlsFingerprint <- function(X, Y, nruns.X, cov, precision = FALSE) {
  X <- as.matrix(X)
  if (! precision) {
    tmpMat <- eigen(cov)
    cov.sinv <- Re(tmpMat$vectors %*% diag(1 / sqrt(tmpMat$values)) %*% 
                     t(tmpMat$vectors))
    ## cov.sinv <- Re(sqrtm(cov))
  } else {
    tmpMat <- eigen(cov)
    cov.sinv <- Re(tmpMat$vectors %*% diag(sqrt(tmpMat$values)) %*% 
                     t(tmpMat$vectors))
  }
  
  output <- tlsLm(cov.sinv %*% X, cov.sinv %*% Y, nruns.X)
  beta.hat <- output$beta.hat
  ci.estim <- output$ci
  sd.estim <- output$sd
  
  beta.hat <- as.vector(beta.hat)
  names(beta.hat) <- colnames(X)
  ## confidence interval
  rownames(ci.estim) <- colnames(X)
  colnames(ci.estim) <- c("Boot lower bound", "Boot upper bound", 
                          "Norm lower bound", "Norm upper bound")
  ## residual consistency test
  result <- list(coefficient = beta.hat,
                 confidence.interval = ci.estim, 
                 var.est = sd.estim)
  result
}


## estimate via Total least square approach
tlsLm <- function(X, Y, nruns.X, conf.level = 0.95) {
  ## input: 
  ##   X: n*k matrix, including k predictors
  ##   Y: n*1 matrix, the observations
  ##   nruns.X: the number of runs used for computing each columns of X
  ## output: 
  ##   a list containing the estimate and confidence interval of the scaling 
  ##   factors estimate.
  ##   coefficient: a k vector, the scaling factors best-estimates
  ##   confidence interval: the lower and upper bounds of the confidence 
  ##   interval on each scaling factor. 
  ##   dcons: the variable used in the residual consistency check. 
  ##   X.tilde: a k*n matrix, the reconstructed responses patterns TLS fit,
  ##   Y.tilde: a 1*n matrix, the reconstructed observations TLS fit.
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
  if(length(nruns.X) != 1) {
    Dn.X <- diag(sqrt(nruns.X))
  } else {
    Dn.X <- sqrt(nruns.X)
  }
  Estls <- function(X, Y, Dn.X) {
    M <- cbind(X, Y)
    lambda <- tail(eigen(t(M) %*% M)$values, 1)
    sigma2.hat <- lambda / n
    Delta.hat <- (t(X) %*% X - lambda * diag(m)) / n
    
    ## beta.hat for the tls regression for adjusted X and Y with equal variance
    beta.hat1 <- as.vector(solve(t(X) %*% X - lambda * diag(m)) %*% t(X) %*% Y)
    ## [I|beta.hat1]
    I.b <- cbind(diag(m), beta.hat1)
    ## var.hat for beta.hat1
    Var.hat1 <- sigma2.hat * (1 + sum(beta.hat1^2)) * 
      (solve(Delta.hat) + sigma2.hat * solve(Delta.hat) %*% 
         solve(I.b %*% t(I.b)) %*% solve(Delta.hat)) / n
    
    ## beta.hat and var.hat for the un prewhitening X and Y
    beta.hat <- beta.hat1 %*% Dn.X
    
    Var.hat <- diag(Var.hat1) * nruns.X
    
    list(beta.hat = beta.hat, Var.hat = Var.hat)
  }
  
  tmp.res <- Estls(X, Y, Dn.X)
  
  beta.hat <- tmp.res$beta.hat
  
  var.hat <- tmp.res$Var.hat
  ## standard deviation from normal approximation
  sd.norm <- sqrt(var.hat)
  ## confidence interval from asymptotical normal distribution
  ci.norm <- cbind(t(beta.hat - 1.96 * sqrt(var.hat)), 
                   t(beta.hat + 1.96 * sqrt(var.hat)))
  colnames(ci.norm) <- c("0.025", "0.975")
  ## nonparamatric bootstrap
  B <- 1000
  for(i in 1:1000) {
    beta.s <- tryCatch({
      resample <- sapply(1:B,
                         function(x) {
                           sample(1:n, size = n, replace = TRUE)
                         })
      beta.s <- apply(resample, 2,
                      function(x) {
                        Xs <- X[x, ]
                        Ys <- Y[x, ]
                        Estls(Xs, Ys, Dn.X)$beta.hat
                      })
      if(ncol(X) == 1) {
        matrix(beta.s, ncol = 1000)
      } else {
        beta.s
      }
    }, error = function(e) {
      as.matrix(c(0, 0))
    })
    if (ncol(beta.s) == B) {
      break
    }
  }
  
  alpha <- 1 - conf.level
  
  ci.ordboot <- t(apply(beta.s, 1, 
                        function(x) {
                          ## x <- x[-which(x > quantile(x, 0.999) | x < quantile(x, 0.001))]
                          quantile(x, c(alpha / 2, alpha / 2 + conf.level))
                        }))
  
  list(beta.hat = beta.hat, ci = cbind(ci.ordboot, ci.norm), sd = sd.norm)
}

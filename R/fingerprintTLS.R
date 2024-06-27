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
##' @param ctlruns a group of independent control runs for estimating 
##' covariance matrix, which is used in two stage bootstrap and the parametric 
##' bootstrap calibration.
##' @param precision indicator for precision matrix, if precision 
##' matrix estimate is used, precision should be set to TRUE.
##' @param conf.level confidence level for confidence interval estimation.
##' @param conf.method method for calibrating the confidence intervals, including
##' no calibration (none), two stage bootstrap (TSB), and parametric bootstrap calibration (PBC).
##' @param cov.method method for estimation of covariance matrix in confidence 
##' interval procedure. It should be consistent to the method to get cov. (only valid if TSB or PBC is considered.)
##' @param B number of replicates in two stage bootstrap and/or parametric bootstrap calibration, 
##' default value is 1000. (only valid if TSB or PBC is considered.)
##' @return a list of the fitted model including point estimate and
##' interval estimate of coefficients and corresponding estimate of 
##' standard error.
##' @author Yan Li
##' @keywords regression tls fingerprinting
##' @references \itemize{ 
##' \item  Gleser (1981), Estimation in a Multivariate "Errors in Variables" 
##' Regression Model: Large Sample Results, \emph{Ann. Stat.} 9(1) 24--44.
##' \item Golub and Laon (1980), An Analysis of the Total Least Squares Problem,
##' \emph{SIAM J. Numer. Anal}. 17(6) 883--893.
##' \item Pesta (2012), Total least squares and bootstrapping with 
##' applications in calibration, \emph{Statistics} 47(5), 966--991.
##' \item Li et al (2021), Uncertainty in Optimal Fingerprinting is Underestimated, \emph{Environ. Res. Lett.} 16(8) 084043.}
##' @noRd

fingerprintTLS <- function(X, Y, cov, nruns.X, ctlruns,
                           precision = FALSE,
                           conf.level = 0.90,
                           conf.method = c("none", "TSB", "PBC", "both"),
                           cov.method = c("l2", "mv"),
                           B = 1000) {
  X <- as.matrix(X)
  cov.method <- match.arg(cov.method)  ## method for the covariance matrix extimation
  conf.method <- match.arg(conf.method)  ## the method for estimating covariance matrix
  if(! conf.method %in% c("none", "TSB", "PBC", "both")) {
    stop("Unknow method for confidence interval construction")
  }
  if(is.null(colnames(X))) {
    colnames(X) <- paste0("forcings ", 1:ncol(X))
  }
  if (! precision) {
    tmpMat <- eigen(cov)
    cov.sinv <- Re(tmpMat$vectors %*% diag(1 / sqrt(tmpMat$values)) %*% 
                     t(tmpMat$vectors))
  } else {
    tmpMat <- eigen(cov)
    cov.sinv <- Re(tmpMat$vectors %*% diag(sqrt(tmpMat$values)) %*% 
                     t(tmpMat$vectors))
  }
  
  output <- tlsLm(cov.sinv %*% X, cov.sinv %*% Y, nruns.X, conf.level)
  beta.hat <- as.vector(output$beta.hat)
  ci.estim <- output$ci
  sd.estim <- output$sd
  
  ## label the cols and rows of the estimation
  numbeta <- length(beta.hat)
  if(is.null(colnames(X))) {
    Xlab <- paste0("forcings ", 1:numbeta)
  } else {
    Xlab <- colnames(X)
  }
  names(beta.hat) <- Xlab
  rownames(ci.estim) <- c(paste0("B: ", Xlab), 
                          paste0("N: ", Xlab))
  colnames(ci.estim) <- paste0(c("Lower ", "Upper "), c((1 - conf.level) / 2, (1 + conf.level) / 2) * 100, "%")
  rownames(sd.estim) <- c("B", "N")
  colnames(sd.estim) <- Xlab
  
  if (conf.method == "TSB" | conf.method == "both") {
    boot <- lapply(1:B,
                   function(i) {
                     for(i in 1:500) {
                       cov.sinv <- tryCatch({
                         resample <- sample(1:nrow(ctlruns),
                                           nrow(ctlruns), replace = TRUE)
                         resample <- unique(resample)
                         ## boot.sample <- genClt(nrow(ctlruns), B = 1, cov)[[1]]
                         boot.sample <- ctlruns[resample, ]
                         if(cov.method == "l2") {
                           cov.boot <- Covest(boot.sample, method = "l2")[[1]]
                         } else if (cov.method == "mv") {
                           cov.boot <- Covest(boot.sample, method = "mv", bandwidth = 0.35)[[1]]
                         }
                         ## compute inverse cov
                         if (! precision) {
                           tmpMat <- eigen(cov.boot)
                           cov.sinv <- Re(tmpMat$vectors %*% diag(1 / sqrt(tmpMat$values)) %*% 
                                            t(tmpMat$vectors))
                         } else {
                           tmpMat <- eigen(cov.boot)
                           cov.sinv <- Re(tmpMat$vectors %*% diag(sqrt(tmpMat$values)) %*% 
                                            t(tmpMat$vectors))
                         }
                         cov.sinv
                       }, error = function(e) {
                         ""
                       })
                       if(cov.sinv[[1]] != "") {
                         break
                       }
                     }
                     output.b <- tlsLm.boot(cov.sinv %*% X, cov.sinv %*% Y, nruns.X, B=1000)
                     output.b
                   })
    
    boot.res <- vector("list", 2)
    
    for(m in 1:length(boot)) {
      boot.res[[1]] <- rbind(boot.res[[1]], boot[[m]]$beta.s)
      boot.res[[2]] <- rbind(boot.res[[2]], t(boot[[m]]$beta.s.list))
    }
    
    FSB.sd <- apply(boot.res[[1]], 2, sd)
    
    ## get critical value
    alpha <- 1 - conf.level
    Z.crt <- qnorm(alpha / 2, lower.tail = FALSE)
    
    ## first stage bootstrap for benchmark
    ci.FSB <- cbind(beta.hat - Z.crt * FSB.sd, 
                    beta.hat + Z.crt * FSB.sd)
    ci.estim <- rbind(ci.estim, ci.FSB)
    sd.estim <- rbind(sd.estim, FSB.sd)
    
    ## two stage bootstrap results
    TSB.sd <- apply(boot.res[[2]], 2, sd)
    ci.TSB <- t(apply(boot.res[[2]], 2, 
                      function(x) {
                        quantile(x, c(alpha / 2, alpha / 2 + conf.level))
                      }))
    ci.estim <- rbind(ci.estim, ci.TSB)
    sd.estim <- rbind(sd.estim, TSB.sd)
    ci.TSBv <- cbind(beta.hat - Z.crt * TSB.sd, 
                     beta.hat + Z.crt * TSB.sd)
    ci.estim <- rbind(ci.estim, ci.TSBv)
    ## confidence interval
    rownames(ci.estim) <- c(paste0("B: ", Xlab), 
                            paste0("N: ", Xlab),
                            paste0("FSB: ", Xlab), 
                            paste0("TSB: ", Xlab), 
                            paste0("TSBv: ", Xlab))
    colnames(sd.estim) <- Xlab
    rownames(sd.estim) <- c("B",
                            "N", 
                            "FSB", 
                            "TSB")
  }
  
  if(conf.method == "PBC" | conf.method == "both") {
    ## compute the ratio for the parametric calibration bootstrap
    
    #### get the fitted X and Y in the model
    output <- tlsFit(cov.sinv %*% X, cov.sinv %*% Y, nruns.X)
    output <- Re(tmpMat$vectors %*% diag(sqrt(tmpMat$values)) %*% 
                   t(tmpMat$vectors)) %*% output
    #### conduct the calibration bootstrap
    #### get number of control runs to estimate the covariance matrix
    rep.num <- nrow(ctlruns)
    #### do the boostrap
    out.beta <- out.var <- NULL
    for (i in 1:B) {
      ## generate new dataset
      for(er in 1:50) {
        tmp.new <- tryCatch({
          Y.new <- MASS::mvrnorm(n = 1, mu = output[, "Y"], Sigma = cov)
          X.new <- NULL
          for(X.ind in 1:ncol(X)) {
            X.new <- cbind(X.new,
                           mvrnorm(n = 1, mu = output[, X.ind], Sigma = cov / nruns.X[X.ind]))
          }
          ctlruns.new <- genClt(rep.num, 1, cov)[[1]]
          if (cov.method == "l2") {
            Cov.new <- Covest(ctlruns.new, method = "l2")$output
          } else if (cov.method == "mv") {
            Cov.new <- Covest(ctlruns.new, method = "mv", bandwidth = 0.35)$output
          }
          ## repeat the function without TSB and PBC
          fingerprintTLS(X.new, Y.new, Cov.new, nruns.X, ctlruns.new, conf.method = "none", cov.method = cov.method)
        }, error = function(e) {
          ""
        })
        if (tmp.new[1] != "") {
          break
        }
      }
      out.beta <- rbind(out.beta, tmp.new$beta)
      out.var <- rbind(out.var, as.vector(tmp.new$var))
    }
    ## return the ratio for the parametric calibration bootstrap
    ratio <- rbind(apply(out.beta, 2, sd), apply(out.beta, 2, sd)) / matrix(colMeans(out.var, na.rm = TRUE), nrow = 2)
  }
  
  ## collect the output
  if(conf.method %in% c("PBC", "both")) {
    result <- list(beta = beta.hat,
                   ci = ci.estim, 
                   var = sd.estim, 
                   pbc.ratio = ratio)
  } else {
    result <- list(beta = beta.hat,
                   ci = ci.estim, 
                   var = sd.estim)
  }
  result
}

## estimate via Total least square approach
tlsLm <- function(X, Y, nruns.X, conf.level) {
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
    ## beta.hat1 <- as.vector(RegGTLS(Y, X))
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
  
  Z.crt <- qnorm((1 - conf.level) / 2, lower.tail = FALSE)
  
  ## confidence interval from asymptotical normal distribution
  ci.norm <- cbind(t(beta.hat - Z.crt * sqrt(var.hat)), 
                   t(beta.hat + Z.crt * sqrt(var.hat)))
  
  colnames(ci.norm) <- c(0.5 - conf.level / 2, 0.5 + conf.level / 2)
  
  B <- 1000
  
  ## nonparamatric bootstrap
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
        matrix(beta.s, ncol = B)
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
                          ## x <- rmout(x)
                          quantile(x, c(alpha / 2, alpha / 2 + conf.level))
                        }))
  
  sd.ordboot <- t(apply(beta.s, 1, 
                        function(x) {
                          ## x <- x[-which(x > quantile(x, 0.99) | x < quantile(x, 0.01))]
                          ## x <- rmout(x)
                          sd(x)
                        }))
  
  sd.norm <- sqrt(var.hat)
  
  ## residual bootstrap
  
  list(beta.hat = beta.hat, ci = rbind(ci.ordboot, ci.norm), 
       sd = rbind(sd.ordboot, sd.norm))
}

## estimate via Total least square approach
## for two stage bootstrap
tlsLm.boot <- function(X, Y, nruns.X, B = 100) {
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
    ## beta.hat1 <- as.vector(RegGTLS(Y, X))
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
  
  ## nonparamatric bootstrap
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
        matrix(beta.s, ncol = B)
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
  
  list(beta.s.list = beta.s, beta.s = Estls(X, Y, Dn.X)$beta.hat)
}

## get the fitted expected response X and Y of the data
tlsFit <- function(X, Y, nruns.X) {
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
    colnames(output) <- c(colnames(X), "Y")
    output
  }
  output <- Estls(X, Y)
  output[, colnames(X)] <- output[, colnames(X)] / 
    t(sqrt(nruns.X) %*% matrix(1, 1, n))
  output
}

## generate independent replicates
genClt <- function(n, B, Cov) {
  lapply(1:B, 
         function(x) {
           MASS::mvrnorm(n, mu = rep(0, nrow(Cov)), Sigma = Cov)
         })
}
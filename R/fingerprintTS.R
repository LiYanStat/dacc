###### internal function for the two sample approach

fingerprintTS <- function(X, Y, nruns.X, cov, Z.2, precision = FALSE, conf.level = 0.90, 
                          B = 1000) {
  ## Z.2 is the second sample for the variance estimation
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
  
  ## compute the estimation using internal function tlsLm
  output <- tlsLmTS(cov.sinv %*% X, cov.sinv %*% Y, nruns.X, conf.level, Z.2, cov.sinv, B)
  beta.hat <- output$beta.hat
  ci.estim <- output$ci
  sd.estim <- output$sd
  
  ## label the cols and rows of the estimation
  numbeta <- length(beta.hat)
  if(is.null(colnames(X))) {
    Xlab <- paste0("forcings ", 1:numbeta)
  } else {
    Xlab <- colnames(X)
  }
  
  rownames(ci.estim) <- c(paste("B", Xlab), paste("N", Xlab))
  colnames(ci.estim) <- c("Lower bound", "Upper bound")
  rownames(sd.estim) <- c("B", "N")
  colnames(sd.estim) <- Xlab
  
  beta.hat <- as.vector(beta.hat)
  names(beta.hat) <- Xlab
  ## combine the results
  result <- list(beta = beta.hat,
                 var = sd.estim,
                 ci = ci.estim)
  return(result)
}

## estimate via Total least square approach
tlsLmTS <- function(X, Y, nruns.X, conf.level, Z.2, cov.sinv, B) {
  ## input: 
  ##   X: n*k matrix, including k predictors
  ##   Y: n*1 matrix, the observations
  ##   nruns.X: the number of runs used for computing each columns of X
  ##   Z.2: the second sample to estimate the variance
  ##   cov.sinv: weight matrix
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
    
    # sigma2.hat <- lambda / n
    # 
    # Delta.hat <- (t(X) %*% X - lambda * diag(m)) / n
    
    ## beta.hat for the tls regression for adjusted X and Y with equal variance
    beta.hat1 <- as.vector(solve(t(X) %*% X - lambda * diag(m)) %*% t(X) %*% Y)
    
    Lambda2 <- diag(t(svd(M)$u) %*% M %*% t(M) %*% svd(M)$u) / diag(t(svd(M)$u) %*% cov.sinv %*% cov(Z.2) %*% cov.sinv %*% svd(M)$u)
    
    XtX <- (svd(M)$v %*% diag(Lambda2) %*% t(svd(M)$v))[1:m, 1:m]
    
    Delta.hat2 <- (XtX - tail(Lambda2, 1) * diag(m)) / n
    
    sigma2.hat2 <- tail(Lambda2, 1) / n
    
    ## [I|beta.hat1]
    I.b <- cbind(diag(m), beta.hat1)
    ## var.hat for beta.hat1
    # Var.hat1 <- sigma2.hat * (1 + sum(beta.hat1^2)) *
    #   solve(Delta.hat) %*% (Delta.hat2 + sigma2.hat *  solve(I.b %*% t(I.b))) %*% solve(Delta.hat) / n
    
    Var.hat1 <- sigma2.hat2 * (1 + sum(beta.hat1^2)) *
      solve(Delta.hat2) %*% (Delta.hat2 + sigma2.hat2 *  solve(I.b %*% t(I.b))) %*% solve(Delta.hat2) / n
    
    ## beta.hat and var.hat for the un prewhitening X and Y
    beta.hat <- beta.hat1 %*% Dn.X
    
    Var.hat <- diag(Var.hat1) * nruns.X
    
    list(beta.hat = beta.hat, Var.hat = Var.hat)
  }
  
  Estls.beta <- function(X, Y, Dn.X) {
    M <- cbind(X, Y)
    lambda <- tail(eigen(t(M) %*% M)$values, 1)
    ## beta.hat for the tls regression for adjusted X and Y with equal variance
    beta.hat1 <- as.vector(solve(t(X) %*% X - lambda * diag(m)) %*% t(X) %*% Y)
    ## beta.hat and var.hat for the un prewhitening X and Y
    beta.hat1 %*% Dn.X
  }
  
  ## compute the estimation
  tmp.res <- Estls(X, Y, Dn.X)
  beta.hat <- tmp.res$beta.hat
  var.hat <- tmp.res$Var.hat
  
  ## get the normal critical value
  Z.crt <- qnorm((1 - conf.level) / 2, lower.tail = FALSE)
  ## confidence interval from asymptotical normal distribution
  ci.norm <- cbind(t(beta.hat - Z.crt * sqrt(var.hat)), 
                   t(beta.hat + Z.crt * sqrt(var.hat)))
  colnames(ci.norm) <- c(0.5 - conf.level / 2, 0.5 + conf.level / 2)
  
  ## compute the normal standard deviation
  sd.norm <- sqrt(var.hat)
  
  ## compute the confidence interval from bootstrap
  ## B <- 1000
  ## nonparamatric bootstrap
  ## error control with tryCatch function
  for(i in 1:10) {
    beta.s <- tryCatch({
      resample <- sapply(1:B,
                         function(x) {
                           sample(1:n, size = n, replace = TRUE)
                         })
      beta.s <- apply(resample, 2,
                      function(x) {
                        Xs <- X[x, ]
                        Ys <- Y[x, ]
                        Estls.beta(Xs, Ys, Dn.X)
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
  
  ## alpha value
  alpha <- 1 - conf.level
  
  ## compute the bootstrap confidence interval
  ci.ordboot <- t(apply(beta.s, 1, 
                        function(x) {
                          x <- rmout(x)
                          quantile(x, c(alpha / 2, alpha / 2 + conf.level))
                        }))
  ## compute the bootstrap standard deviation
  sd.ordboot <- t(apply(beta.s, 1, 
                        function(x) {
                          x <- rmout(x)
                          sd(x)
                        }))
  
  list(beta.hat = beta.hat, ci = rbind(ci.ordboot, ci.norm), 
       sd = rbind(sd.ordboot, sd.norm))
}

## functions for removing possible outliers
rmout <- function(x) {
  ind <- c(which(x %in% boxplot.stats(x)$out),  which(is.na(x)))
  ind <- unique(ind)
  if(length(ind) == 0) {
    x
  } else {
    x[-ind]
  }
}

rmoutInt <- function(x) {
  ind <- c(which(x[, 1] %in% boxplot.stats(x[, 1])$out), 
           which(x[, 2] %in% boxplot.stats(x[, 2])$out), 
           which(is.na(x[, 1])), 
           which(is.na(x[, 2])))
  ind <- unique(ind)
  if(length(ind) == 0) {
    x
  } else {
    x[-ind, , drop = FALSE]
  }
}
##' Regularized Optimal Fingerprinting via linear regression.
##'
##' This function detects the signal factors on the observed data via linear 
##' regression. 
##'
##' @param Y observed data.
##' @param X signal pattern to be detected.
##' @param nruns.X the number of runs used for computing each columns of X
##' @param Proj projection matrix indicating the relationship of the signal 
##' used in regression and the signal to be reported.
##' @param signal names of the external forces to be detected.
##' @param lm.method the Least sqaure method to be used, can be "OLS" or "TLS".
##' @param Z.pw a matrix containing sample from contorl runs to estimate 
##' the covariance matrix in prewhitening.
##' @param Z.rct amatrix containing sample from contorl runs to estimate the 
##' covariance matrix in residual consistency checking.
##' @param rct.method test used to check the residual consistence, 
##' method should be "OLS_AT99" or "OLS_Corr" for "OLS" lm.method 
##' and "MC" or "AS03" for "TLS".
##' @param confidence.level confidence level for OLS method.
##' @param conf.level confidence level for TLS method.
##' @param ci.approach formula for computing confidence interval in TLS.
##' @param cov.method algorithm for computing the covariance matrix.
##' @return a list of the fitted model including point estimate and 
##' interval estimate of coefficient, p value of the residual consistency test.
##' @author Yan Li
##' @references Ribes et al. 2013, Allen and Stott 2003, Fan et al. 2015, 
##' Friedman et al 2008.
##' @keywords Regularized Fingerprinting, OLS, TLS
##' @examples
##' data(cnrmDat)
##' Y <- cnrmDat$obs
##' st.n <- length(Y)
##' nZ <- length(cnrmDat$ctlruns) / st.n
##' Z <- matrix(cnrmDat$ctlruns, nrow = st.n, ncol = nZ) 
##' ## extract Z
##' ex.Z <- extractZ(Z, frac.Z = 0.5)
##' Z.pw <- ex.Z$Z1
##' Z.rct <- ex.Z$Z2
##' ## TLS for 2 factors case
##' factors <- c("ANT", "NAT")
##' ## Projection matrix
##' Proj <- diag(length(factors))
##' X <- cnrmDat$signal[, factors]
##' X <- as.matrix(X)
##' nruns.X <- c(10, 6)
##' ## different residual consistency test 
##' regFingerprint(Y, X, nruns.X, Proj, factors, lm.method = "TLS", Z.pw, Z.rct, 
##' rct.method = "MC")
##' @import MASS stats utils glasso flare methods
##' @importFrom expm sqrtm
##' @export regFingerprint
regFingerprint <- function(Y, X, nruns.X, Proj, signal, lm.method, Z.pw, Z.rct, rct.method, 
                           confidence.level = 0.05, 
                           conf.level = 0.1, 
                           ci.approach = c("AS03", "ODP"), 
                           cov.method = c("LW", "glasso", "flare", "fastclime")) {
  X <- as.matrix(X)
  cov.method <- match.arg(cov.method)
  ## compute the covariance matrix and precision matrix
  if (cov.method == "LW") {
    cov <- lwRegcov(t(Z.pw))
    cov.inv <- solve(cov)
  } else if (cov.method == "glasso") {
    S <- var(t(Z.pw))
    tmpList <- glasso(S, rho=.02)
    cov <- tmpList$w
    cov.inv <- tmpList$wi
  } else if (cov.method == "flare") {
    out <- sugm(t(Z.pw), verbose = FALSE)
    output <- sugm.select(out, verbose = FALSE)
    cov.inv <- as.matrix(output$opt.icov)
    cov <- solve(cov.inv)
  } else if (cov.method == "fastclime") {
    cov <- NULL
    cov.inv <- NULL
  } else {
    stop("unknown algorithm for estimating covariance and precision matrix")
  }
  ## estimate coefficient and confidence interval
  if (lm.method == "OLS") {
    ## l <- lm.gls(Y ~ X - 1, W = solve(cov))
    Ft <- t(solve(t(X) %*% cov.inv %*% X) %*% t(X) %*% cov.inv)
    beta.hat <- t(Y) %*% Ft %*% t(Proj)  ## estimate of coefficients
    ## compute the covariance matrix and confidence interval of estimate
    Cov.valid = Z.rct %*% t(Z.rct) / ncol(Z.rct)
    ## covariance matrix
    cov.betahat <- Proj %*% t(Ft) %*% Cov.valid %*% Ft %*% t(Proj)
    ## confidence interval
    betaHat.inf <- beta.hat - 
      qt(1 - confidence.level / 2, ncol(Z.rct)) * sqrt(t(diag(cov.betahat)))
    betaHat.sup <- beta.hat + 
      qt(1 - confidence.level / 2, ncol(Z.rct)) * sqrt(t(diag(cov.betahat)))
    ## residual consistency checking
    ## compute the pseudo inverse of the 
    ginv.valid <- ginv(Cov.valid)
    df <- (length(Y) - ncol(X))
    ## compute the residual
    epsilon <- Y - X %*% solve(Proj) %*% t(beta.hat)
    if (rct.method == "OLS_AT99") {  ## Formula provided by Allen & Tett (1999)
      dcons <- t(epsilon) %*% ginv.valid %*% epsilon / df
      pvalue <- 1 - pf(dcons, df1 = df, df2 = ncol(Z.rct))
    } else if (rct.method == "OLS_Corr") {
      df <- (length(Y) - ncol(X))
      if (ncol(Z.rct) - length(Y) + 1 > 0) {
        dcons <- t(epsilon) %*% ginv.valid %*% epsilon / (ncol(Z.rct) * df) * 
          (ncol(Z.rct) - length(Y) + 1)
        pvalue <- 1 - pf(dcons, df1 = df, df2 = ncol(Z.rct) - length(Y) + 1)
      } else {
        pvalue <- NaN
      }
    } else {
      pvalue <- "Unknown Consistency test"
    }
    ## format the output
    ## estimates of coefficient
    beta.hat <- as.vector(beta.hat)
    names(beta.hat) <- signal
    ## covariance matrix
    rownames(cov.betahat) <- colnames(cov.betahat) <- signal
    ## confidence interval
    ci.estim <- t(rbind(betaHat.inf, betaHat.sup))
    rownames(ci.estim) <- signal
    colnames(ci.estim) <- c("lower bound", "upper bound")
    ## residual consistency test
    consistency.test <- pvalue
    consistency.test <- as.vector(consistency.test)
    names(consistency.test) <- "p value"
    result <- list(coefficient = beta.hat,
                   covariance = cov.betahat,
                   confidence.interval = ci.estim,
                   consistency.test = consistency.test)
  } else if (lm.method == "TLS") {
    cov.sinv <- sqrtm(cov.inv)
    cov.s <- sqrtm(cov)
    ci.approach <- match.arg(ci.approach)
    output <- tlsLm(t(X) %*% cov.sinv, t(Y) %*% cov.sinv, t(Z.rct)  %*% cov.sinv, 
                    nruns.X, Proj, ci.approach)
    ## reconstructed data set
    X.tilde <- cov.s %*% t(output$X.tilde)
    Y.tilde <- cov.s %*% output$Y.tilde
    ## estimate and confidence interval of coefficient
    beta.hat <- output$beta.hat
    dcons <- output$dcons
    betaHat.inf <- output$betaHat.inf
    betaHat.sup <- output$betaHat.sup
    ## residual consistency check
    n.pw <- dim(Z.pw)[2]
    n.rct <- dim(Z.rct)[2]
    if (rct.method == "MC") {  ## Monte-Carlo simulations
      dcon.h0 <- conMcTls(cov, X, beta.hat, nruns.X, n.pw, n.rct, ci.approach)
      ## gaussian kernal estimate
      B <- length(dcon.h0)
      h <- 1.06 * sd(dcon.h0) * B^(-1 / 2)
      pvalue <- mean(1 - pnorm(dcons * rep(1, B), dcon.h0, h * rep(1, B)))
    } else if (rct.method == "AS03") {
      pvalue <- 1 - pf(dcons / (length(Y) - ncol(X)), (length(Y) - ncol(X)), n.rct)
    } else {
      stop("Unknow approach for residual consistency test")
    }
    ## format the output
    ## estimates of coefficient
    beta.hat <- as.vector(beta.hat)
    names(beta.hat) <- signal
    ## confidence interval
    ci.estim <- t(rbind(betaHat.inf, betaHat.sup))
    rownames(ci.estim) <- signal
    colnames(ci.estim) <- c("lower bound", "upper bound")
    ## residual consistency test
    consistency.test <- pvalue
    consistency.test <- as.vector(consistency.test)
    names(consistency.test) <- "p value"
    result <- list(coefficient = beta.hat,
                   confidence.interval = ci.estim,
                   consistency.test = consistency.test)
  } else {
    stop("approach for fitting linear model should be 'OLS' or 'TLS'")
  }
  result
}

## estimate via Total least square approach
tlsLm <- function(X, Y, Z, nruns.X, Proj, ci.approach, 
                  conf.level = 0.1) {
  ## input: 
  ##   X: k*n matrix, including k predictors
  ##   Y: 1*n matrix, the observations
  ##   Z: a set of independant realisations of internal variability.
  ##   nruns.X: the number of runs used for computing each columns of X
  ##   Proj: the matrix that defined which forcings are used in each simulations.
  ##   ci.approach: a string, the formula that will be used to compute the 
  ##   confidence intervals; "AS03" is recommended, the alternative is "ODP".
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
    stop("Y should be a 1*n matrix")
  }
  if (dim(X)[2] != dim(Y)[2])  {  ## check size of X and Y
    stop("sizes of inputs X, Y are not consistent")
  }
  n <- dim(Y)[2]  ## number of observations
  m <- dim(X)[1]  ## number of predictors
  ## Normalise the variance of X
  X <- X * sqrt(nruns.X) %*% matrix(1, 1, n)
  Dn.X <- diag(sqrt(nruns.X))
  ## compute the estimate of beta via svd
  M <- rbind(X, Y)
  ## svd decomposition
  U <- svd(M)$u
  D <- svd(M)$d
  V <- svd(M)$v
  ## get the last singular vector
  U.last <- U[, ncol(U)]
  U.proj <- rbind(Proj %*% Dn.X %*% U.last[1:(length(U.last)-1)], 
                  U.last[length(U.last)])
  ## compute the beta_hat
  beta.hat <- - U.proj[-length(U.proj)] / U.proj[length(U.proj)]
  ## reconstruct data
  D.tilde <- diag(D)
  D.tilde[m + 1, m + 1] <- 0
  Z.tilde <- U %*% D.tilde %*% t(V)
  X.tilde <- Z.tilde[1:m, ] / sqrt(nruns.X) %*% matrix(1, 1, n)	
  Y.tilde <- Z.tilde[m+1, ]
  
  ## compute the confidence interval
  d <- D^2  ## The square singular values (denoted by lambda in AS03)
  ## compute corrected singular value (cf Eq 34 in AS03)
  nZ <- dim(Z)[1]
  d.hat <- NULL
  for (i in 1:length(d)) {
    vi <- t(V[, i])
    if (ci.approach == "AS03") {
      d.hat <- c(d.hat, d[i] / ((vi %*% t(Z) %*% Z / nZ) %*% t(vi)))
    } else if (ci.approach == "ODP") {
      d.hat <- c(d.hat, d[i] / (vi^2 %*%  colSums(Z^2) / nZ))
    } else {
      stop("unknown formula for computation of TLS CI")
    }
  }
  ## the "last" corrected singular value will be used in the 
  ## Residual Consistency Check
  dcons <- tail(d.hat, 1)
  ## compute the threshold of the F distribution for CI estimation
  thres.f <- sqrt(qf(1 - conf.level, df1 = 1, df2 = nZ))
  ## in order to compute CI, we need to run through the m-sphere (cf Eq 30 in AS03)
  npt <- 1000  ## number of points on the sphere
  if (m == 1) {
    Pts <- rbind(1, -1)
  } else {
    Pts.R <- sapply(1:m, 
                    function(x) {
                      rnorm(npt, 0, 1)
                    })
    Pts <- Pts.R / (sqrt(rowSums(Pts.R^2)) %*% matrix(1, 1, m))
  }
  deltaD.hat <- d.hat - min(d.hat)
  a <- thres.f * Pts
  if (sum(head(deltaD.hat, length(deltaD.hat))) != 0) {
    b.m1 <- a / (matrix(1, nrow = nrow(Pts), ncol = 1) %*% 
                   t(sqrt(deltaD.hat[1:(length(deltaD.hat) - 1)])))
    b.m2 <- sqrt(1 - rowSums(b.m1^2))
  } else {  ## the confidence interval will be unbounded
    betaHat.inf = rep(NaN, length(beta.hat))
    betaHat.sup = rep(NaN, length(beta.hat))
  }
  ## b.m2 need to be strctly positive, otherwise the CI will be unbounded
  if (! all(Re(b.m2) == b.m2) | sum(b.m2) == 0) {
    betaHat.inf = rep(NaN, length(beta.hat))
    betaHat.sup = rep(NaN, length(beta.hat))
  } else {
    ## Then in order to CI that include +/- infinity, 
    ## the computation are made in terms of angles, 
    ## based on complex numbers (this is a descrepancy with ODP)
    V.pts <- cbind(b.m1, b.m2) %*% t(U)
    nInd <- ncol(V.pts)
    ## project on the scaling factors to be detected
    Vpts.proj <- cbind(V.pts[, 1:(nInd - 1)] %*% Dn.X %*% t(Proj), V.pts[, nInd])
    betaHat.inf <- betaHat.sup <- NULL
    tmpVec <- t(apply(Vpts.proj, 1, function(x) {- x[-nInd] / x[nInd]}))
    CI <- matrix(NA, nrow = m, ncol = 2)
    for(i in 1:m) {
      betaHat.inf <- c(betaHat.inf, min(tmpVec[, i]))
      betaHat.sup <- c(betaHat.sup, max(tmpVec[, i]))
    }
    ## for (i in 1:m) {
    ##   Vc.2d.pts <- complex(real = Vpts.proj[, i], imaginary = Vpts.proj[, nInd])
    ##   Vc.2d.ref <- complex(real = U.proj[i], imaginary = tail(U.proj, 1))
    ##   Vprod.2d  <- Vc.2d.pts / Vc.2d.ref
    ##   arg <- sort(Im(log(Vprod.2d)))
    ##   deltaArg.min <- head(arg, 1)
    ##   deltaArg.max <- tail(arg, 1)
    ##   delta.max.1 <- max(arg[2:length(arg)] - arg[1:(length(arg) - 1)])
    ##   k1 <- which(arg[2:length(arg)] - arg[1:(length(arg) - 1)] == delta.max.1)
    ##   tmpL <- c(delta.max.1, deltaArg.max - deltaArg.min + 2 * pi)
    ##   delta.max <- max(tmpL)
    ##   k2 <- which(tmpL == delta.max)
    ##   if (delta.max < pi) {
    ##     betaHat.inf = rep(NaN, length(beta.hat))
    ##     betaHat.sup = rep(NaN, length(beta.hat))
    ##   } else {
    ##     if (k2 != 2) {
    ##       warning("warning k2")
    ##     }
    ##     arg.ref <- Im(log(Vc.2d.ref))
    ##     arg.min <- deltaArg.min + arg.ref
    ##     arg.max <- deltaArg.max + arg.ref
    ##     betaHat.inf <- c(betaHat.inf, -1 / tan(arg.min))
    ##     betaHat.sup <- c(betaHat.sup, -1 / tan(arg.max))
    ##   }
    ## }
  }
  list(beta.hat = beta.hat, 
       betaHat.inf = betaHat.inf, 
       betaHat.sup = betaHat.sup, 
       dcons = dcons, 
       X.tilde = X.tilde, 
       Y.tilde = Y.tilde)
}

## Monte-Carlo simulations to get the empirical null distribution of the 
## lambda for residual consistency check
conMcTls <- function(cov, X, beta, nruns.X, n.pw, n.rct, ci.approach, MC = 1000) {
  ## input:
  ##   cov: covariance matrix of the observation
  ##   X: a n*I matrix, the "true" predictors used in the regression equation.
  ##   nruns.X: a I vector, the number of runs used for evaluating 
  ##   each response pattern X.
  ##   n.pw: sample size of the replicates used in prewhitening.
  ##   n.rct: sample size of the replicates used in the 
  ##   residual consistency statistic computing.
  ##   ci.approach: he formula to be used for computing the RCC variable.
  ##   MC: replication of bootstrap
  if (nrow(cov) != ncol(cov)) {
    ## check the covariance matrix is square matrix
    stop("covariance matrix should be square matrix")
  } 
  n.pred <- ncol(X)  ## number of predictors
  n <- nrow(cov)
  beta0 <- as.vector(beta)
  ## beta0 <- rep(1, n.pred)
  dcons.h0 <- NULL
  for (i in 1:MC) {
    Y.it <- X %*% beta0
    Y.i <- Y.it + mvrnorm(n = 1, mu = rep(0, n), Sigma = cov)
    X.i <- X + sapply(nruns.X, 
                      function(m) {
                        mvrnorm(n = 1, mu = rep(0, n), Sigma = cov / m)
                      })
    X.ic <- rep(1, n) %*% t((sqrt(nruns.X))) * X.i
    Z1 <- t(mvrnorm(n = n.pw, mu = rep(0, n), Sigma = cov))
    Z2 <- t(mvrnorm(n = n.rct, mu = rep(0, n), Sigma = cov))
    ## compute the covariance
    cov.hat <- lwRegcov(t(Z1))
    cov.shat <- sqrtm(solve(cov.hat))
    M <- cov.shat %*% cbind(X.ic, Y.i)
    U.i <- svd(t(M))$u
    V.i <- svd(t(M))$v
    D.i <- svd(t(M))$d
    d <- D.i^2
    nd <- length(d)
    Vi <- t(V.i[, nd])
    Z2w <- t(cov.shat %*% Z2)
    if (ci.approach == "AS03") {
      dcons.h0 <- c(dcons.h0, d[nd] / (Vi %*% (t(Z2w) %*% Z2w / n.rct) %*% t(Vi)))
    } else if (ci.approach == "ODP") {
      dcons.h0 <- c(dcons.h0, d[nd] / (Vi^2 %*% colSums(Z2w^2) / n.rct))
    }
  }
  dcons.h0
}

##' Pre-processing of data set via weighting of spherical harmonics
##'
##' This function processes the raw data via spherical harmonics and full 
##' rank projection.
##' 
##' @param X matrix containing the signal pattern, 
##' @param Y observed data
##' @param Z n*1 matrix containing the output from the control runs
##' @param frac.Z the propotion used to devide the replicates Z into two sets
##' @param trunc Number of truncation used in deviding the spatial areas. 
##' The spatial dimension of X, Y and Z should be (trunc + 1)^2
##' And the total length of Y should be the spatial * temporal dimension
##' Y observed data.
##' @return a list of the processed data sets
##' @author Yan Li
##' @references Ribes et al. 2013.
##' @keywords spherical harmonics
##' @examples
##' data(cnrmDat)
##' Y <- cnrmDat$obs
##' factors <- c("ANT", "NAT")
##' X <- cnrmDat$signal[, factors]
##' X <- as.matrix(X)
##' output <- preProcess(X, Y, cnrmDat$ctlruns, frac.Z = 0.5, trunc = 0)
##' @export preProcess
preProcess <- function(X, Y, Z, frac.Z, trunc) {
  st.n <- length(Y)
  ## check whether the data sets are correctly recorded
  if (st.n != nrow(X)) {
    stop("external forces should be consistent to observed data")
  }
  spa.n <- (trunc + 1)^2
  temp.n <- st.n / spa.n  ## compute the number of time steps
  if (temp.n != floor(temp.n)) {
    stop("error reading the data sets")
  }
  nZ <- length(Z) / st.n
  Z <- matrix(Z, nrow = st.n, ncol = nZ)
  ## extract total sample Z
  ex.Z <- extractZ(Z, frac.Z)
  Z1 <- ex.Z$Z1
  Z2 <- ex.Z$Z2
  ## compute the total wave number under the given resolution of truncation
  twave <- twNumber(trunc)
  p <- 1 / sqrt(2 * twave - 1)
  Pml = diag(as.vector((p %*% t(rep(1, temp.n)))))
  ## convert the original data by weighting process
  y <- Pml %*% Y
  Z1 <- Pml %*% Z1
  Z2 <- Pml %*% Z2
  X <- Pml %*% X
  ## remove useless dimension to make the covariance matrix full rank
  ## equivalent to remove one time step
  red.n <- st.n - spa.n
  U <- projFullank(temp.n, spa.n) 
  ## project all input data
  Y.c <- U %*% Y
  Z1.c <- U %*% Z1
  Z2.c <- U %*% Z2
  X.c <- U %*% X
  list(Y = Y.c, Z1 = Z1.c, Z2 = Z2.c, X = X.c)
}

## Generating total wave number of the spherical harmonics
twNumber <- function(n) {
  ## input: 
  ##   n: Number of truncation, spatial demension is equal to (Trunc + 1)^2
  ## output:
  ##   l: (n+1)^2*1 matrix containing the total wave number.
  nr <- (n + 1)^2
  l <- matrix(0, nrow = nr, ncol = 1)
  l[1:(n+1)] = 1:(n+1)
  ir <- n + 2
  if (n > 0) {
    for (i in 2:(n+1)) {
      for (j in i:(n+1)) {
        l[c(ir,ir+1)] = c(j, j)
        ir <- ir + 2
      }
    }
  }
  l
}

##' This function devides the total replicates Z into two separate sets Z1 and Z2, 
##' where Z1 for prewhiting and Z2 is for residual consistency test.
##' @param Z matrix containing the output from control runs for internal 
##' climate variability.
##' @param frac.Z the propotion used to devide the replicates Z into two sets
##' @param approach approach used to devide the total sample Z. 
##' 'regular', one vector out of two; 'random'  purely random; 
##' 'segment' the first n2 vectors in Z is asigned for Z2.
##' @return a list of the processed data sets Z1 and Z2
##' @author Yan Li
##' @references Ribes et al. 2013.
##' @examples
##' data(cnrmDat)
##' Y <- cnrmDat$obs
##' st.n <- length(Y)
##' nZ <- length(cnrmDat$ctlruns) / st.n
##' Z <- matrix(cnrmDat$ctlruns, nrow = st.n, ncol = nZ)
##' ## extract Z
##' output <- extractZ(Z, frac.Z = 0.5)
##' @export extractZ
extractZ <- function(Z, frac.Z, 
                     approach = c("regular", "random", "segment")) {
  approach <- match.arg(approach)
  if (! approach %in% c("regular", "random", "segment")) {
    stop("apprach should be one of 'regular', 'random' and 'segment'")
  }
  nZ <- ncol(Z)
  nZ.2 <- floor(nZ * frac.Z)
  if (approach == "segment") {
    Z2 <- Z[, 1:nZ.2]
    Z1 <- Z[, -(1:nZ.2)]
  } else if (approach == "random") {
    u <- rnorm(nZ, 0, 1)
    z <- order(u, decreasing = TRUE)
    Z2 <- Z[, z[1:nZ.2]]
    Z1 <- Z[, z[-(1:nZ.2)]]
  } else {
    ind <- floor(seq(1/frac.Z, nZ, by = (1/frac.Z)))
    Z2 <- Z[, ind]
    Z1 <- Z[, -ind]
  }
  list(Z1 = Z1, Z2 = Z2)
}

## compute the eigenvalues and corresponding eigenvectors
speco <- function(M) {
  eigenVa.M <- diag(eigen(M)$values)
  eigenVe.M <- eigen(M)$vectors
  if (max(Im(eigenVa.M)) / max(Re(eigenVa.M))  > 10^(-12)) {
    stop("matrix is not symmetric")
  }
  P1 <- Re(eigenVe.M)
  D1 <- Re(diag(eigenVa.M))
  ind <- order(D1, decreasing = TRUE)
  P <- P1[, ind]
  D <- diag(D1[ind])
  list(P = P, D = D)
}

## generate projection matrix used in transform data to full rank data set
## References : 
## Ribes A., S. Planton, L. Terray (2012) Regularised optimal 
## fingerprint for attribution. Part I: method, 
## properties and idealised analysis. Climate dynamics
projFullank <- function(t.n, spa.n) {
  ## input:
  ##   t.n: the number of time steps.
  ##   spa.n: the number of spatial steps.
  ## output:  
  ##   a (t.n - 1) * spa.n x t.n * spa.n matrix, the projection matrix.
  ## the matrix corresponding to the temporal centering
  M = diag(t.n) - matrix(1, t.n, t.n) / t.n
  eigen.vector <- speco(M)$P
  eigen.value <- speco(M)$D
  Eigen.vector <- eigen.vector[, 1:(t.n - 1)]
  Proj <- matrix(0, (t.n - 1) * spa.n, t.n * spa.n)
  for (i in 1:spa.n) {
    rowInd <- seq(i, ((t.n - 1) * spa.n), by = spa.n)
    colInd <- seq(i, (t.n * spa.n), by = spa.n)
    Proj[rowInd, colInd] <- Eigen.vector
  }
  Proj
}



##' Large-Dimensional Covariance Matrix Estimators
##' This function provide several estimators for the 
##' covaraince matrix under different loss functions
##'
##' @param X sample used to estimate covariance matrix.
##' X should be a n*p matrix indicating sample size n 
##' and dimension of covariance matrix p
##' @param loss loss function to minimize, should be 
##' one of "l2", "stein", "mv"
##' @return a estimate of covariance matrix.
##' @author Yan Li
##' @references 
##' Ledoit O., M. Wolf (2004) A well-conditioned estimator for 
##' large-dimensional covariance matrices. 
##' Journal of Multivariate Analysis, 88(2): 365-411.
##' 
##' Ledoit O., M. Wolf (2017) 
##' Direct Nonlinear Shrinkage Estimation of 
##' Large-Dimensional Covariance Matrice.
##' Working paper. avaliable at 
##' https://ssrn.com/abstract=3047302
##' 
##' Oliver Ledoit and Michael Wolf (2017)
##' Optimal Estimation of a Large-Dimensional 
##' Covariance Matrix under Stein's Loss
##' Working paper. avaliable at 
##' http://dx.doi.org/10.2139/ssrn.2264903
##' 
##' @keywords large-dimension, covariance matrix, stein's loss
##' @examples
##' X <- sapply(1:20, function(x) {rnorm(10, )})
##' Covest(X, loss = "mv")
##' @import Iso
##' @export Covest
Covest <- function(X, loss = c("l2", "stein", "mv")) {
  loss <- match.arg(loss)
  if (!loss %in% c("l2", "stein", "mv")) {
    stop("please select one objective function from 'l2', 'stein' or 'mv'")
  }
  output <- NULL
  if (loss == "l2") {
    ## Linear shrinkage estimator under frobenious loss function
    output <- lwRegcov(X)
  } else if (loss == "mv") {
    ## Nonlinear shrinkage estimator under minimum variance loss function
    output <- lwMvcov(X)
  } else {
    ## Nonlinear shrinkage estimator under stein's loss function
    ## output <- lwSteincov(X)
  }
  output
}

## compute the regularized estimate of the covariance matrix
## References:
## Ledoit O., M. Wolf (2004) A well-conditioned estimator for 
## large-dimensional covariance matrices. 
## Journal of Multivariate Analysis, 88(2): 365-411.
## L&W estimate

lwRegcov <- function(X) {
  ## input: 
  ##   the independent replicates of the random vector, denoted by matrix X
  ##   each column represents a random variable. Thus X is n*p matrix
  ## output:
  ##   regularized covariance matrix
  n <- nrow(X)
  p <- ncol(X)
  sample.cov <- t(X) %*% X / n
  Ip <- diag(p)
  m <- sum(diag(sample.cov %*% Ip)) / p ## first estimate in L&W
  Xp <- sample.cov - m * Ip 
  d2 <- sum(diag(Xp %*% t(Xp))) / p  ## second estimate in L&W
  bt <- NULL
  for (i in 1:n) {
    Xr.i <- matrix(X[i, ], nrow = 1)
    M.i <- t(Xr.i) %*% Xr.i
    bt <- c(bt, sum(diag((M.i - sample.cov) %*% t(M.i - sample.cov))) / p)
  }
  bb2 <- 1 / n^2 * sum(bt)
  b2 <- min(bb2,d2)  ## third estimate in L&W
  a2 <- d2 - b2  ## fourth estimate in L&W
  ## the regularized estimate of the covariance matrix
  b2 * m / d2 * Ip + a2 / d2 * sample.cov
}

## compute the regularized estimate of the covariance matrix
## References:
## Ledoit O., M. Wolf (2017) 
## Direct Nonlinear Shrinkage Estimation of 
## Large-Dimensional Covariance Matrice.
## Working paper. avaliable at 
## https://ssrn.com/abstract=3047302

lwMvcov <- function(X) {
  ## input: 
  ##   the sample replicates of the random vector, 
  ##   denoted by n*p matrix X.
  ## output:
  ##   nonlinear shrinkage covariance matrix under minimum variance 
  ##   loss function.
  n <- nrow(X)
  p <- ncol(X)
  ## sample covariance matrix
  sample.cov <- t(X) %*% X / n
  ## eigen value and vectors
  lambda <- eigen(sample.cov)$value
  eigen.v <- eigen(sample.cov)$vectors
  temOd <- order(lambda)
  lambda <- lambda[temOd]
  eigen.v <- eigen.v[, temOd]
  ## compute the direct kernal estimator
  lambda.tilde <- lambda[max(1, p - n  + 1):p]
  h <- n^-0.35
  L <- matrix(rep(lambda.tilde, min(p, n)), length(lambda.tilde), min(p,n))
  TemM <- 4 * t(L)^2 * h^2 - (L - t(L))^2
  TemM[TemM < 0] = 0
  ftilde <- colMeans(sqrt(TemM) / (2 * pi * h^2 * t(L)))
  TemM <- -(4 * t(L)^2 * h^2 - (L - t(L))^2)
  TemM[TemM < 0] = 0
  Hftilde <- colMeans((sign(L - t(L)) * sqrt(TemM) - L + t(L)) / 
                        (2 * pi * h^2 * t(L)))
  if(p <= n) {
    ## non-singular case
    dtilde <- lambda.tilde / ((pi * p / n * lambda.tilde * ftilde)^2 + 
                                (1 - p/n - pi * p/n * lambda.tilde * Hftilde)^2)
  } else {
    ## singular case
    Hftilde0 <- (1 - sqrt(1 - 4 * h^2)) / (2 * pi * h^2) * mean(1 / lambda.tilde)
    dtilde0 <- 1 / (pi * (p - n) / n * Hftilde0)
    dtilde1 <- lambda.tilde / (pi^2 * lambda.tilde^2 * (ftilde^2 + Hftilde^2))
    dtilde <- c(rep(dtilde0, p-n), dtilde1)
  }
  ## Pool Adjacent Violators algorithm.
  dtilde <- pava(dtilde)
  ## output the covariance matrix
  eigen.v %*% diag(dtilde) %*% t(eigen.v)
}

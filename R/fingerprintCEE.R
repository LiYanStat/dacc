##' Optimal Fingerprinting with Estimating Equations.
##'
##' This function estimates the signal factors and corresponding confidence 
##' interval via the estimating equation methods.
##'
##' @param Xtilde \eqn{n \times p} matrix, signal pattern to be detected.
##' @param Y \eqn{n \times 1} matrix, length \eqn{nS \times nT}, observed data.
##' @param mruns number of ensembles to estimate the corresponding pattern. 
##' It is used as the scale of the covariance matrix for Xi.
##' @param ctlruns.1 \eqn{m \times n} matrix, a group of \eqn{m} independent control 
##' runs for estimating covariance matrix, which is used in point estimation of 
##' the signal factors.
##' @param ctlruns.2 \eqn{m \times n} matrix, another group of \eqn{m} independent control 
##' runs for estimating the corresponding confidence interval of the signal factors, 
##' default same as ctlruns.1.
##' @param nS number of locations for the observed responses.
##' @param nT number of time steps for the observed responses.
##' @param nB number of replicates in bootstrap for the bootstrapped variance and confidence 
##' interval estimations.
##' @param conf.level confidence level for confidence interval estimation.
##' @param cal.a indicator for calculating the a value, otherwise use default value a = 1.
##' @param missing indicator for whether missing values present in Y.
##' @param ridge shrinkage value for adjusting the method for missing observations if missing = TRUE.
##' @return a list of the fitted model including point estimate and
##' interval estimate of beta and a coefficients and corresponding estimate of 
##' variance, together with the residuals for model diagnostics.
##' @author Yan Li
##' @references \itemize{ 
##' \item Sai et al (2023), Optimal Fingerprinting with Estimating Equations,
##'       \emph{Journal of Climate} 36(20), 7109–-7122.}
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
##'                 MASS::mvrnorm(n = 1, mu = NAT, Sigma = Cov / mruns[2]))
##' ## control runs
##' ctlruns <- MASS::mvrnorm(100, mu = rep(0, nrow(Cov)), Sigma = Cov)
##' ## ctlruns1 for the point estimation and ctlruns2 for the interval estimation
##' ctlrun.1 <- ctlrun.2 <- ctlruns
##' ## number of locations
##' nS <- 25
##' ## number of year steps
##' nT <- 10
##' ## call the function to estimate the signal factors via CEE
##' fingerprintCEE(Xtilde, Y, mruns, 
##'                ctlrun.1, ctlrun.2,
##'                nS, nT, nB = 10,
##'                conf.level = 0.9, 
##'                cal.a = TRUE,
##'                missing = FALSE, ridge = 0)
##' @import magrittr
##' @importFrom MASS mvrnorm
##' @importFrom stats cov qnorm quantile sd var pnorm
##' @importFrom utils tail
##' @importFrom pracma ceil
##' @importFrom janitor remove_empty
##' @export fingerprintCEE

fingerprintCEE <- function(Xtilde, Y, mruns,
                           ctlruns.1, ctlruns.2,
                           nS, nT,
                           nB = 0, 
                           conf.level = 0.9,
                           cal.a = TRUE,
                           missing = FALSE,
                           ridge = 0) {
  ## check whether there is missing value in the response
  output <- NULL
  if(missing) {
    ## if missing, do not calculate the a value
    output <- 
      eefp_mis(Xt = Xtilde, Y = Y, m = mruns,
               ctlruns1 = ctlruns.1, ctlruns2 = ctlruns.2,
               ni = nS, C = nT,
               ridge = ridge, nB = nB, conf.level = conf.level)
  } else {
    output <- 
      eefp(Xt = Xtilde, Y = Y, m = mruns,
           ctlruns1 = ctlruns.1, ctlruns2 = ctlruns.2,
           ni = nS, C = nT,
           nB = nB, conf.level = conf.level, cal_a = cal.a)
  }
  return(output)
}

## fingerprinting model by EE method
#### Xt:         n \times p matrix, the ensemble model simulations, each column represents a
####             external forcing.
#### Y:          n \times 1 matrix, length ni (number of locations) x C (number of time steps), 
####             the corresponding observed values.
#### m:          p vector, the number of ensemble runs to obtain the signals of external forcings.
#### ctlruns1:   control runs used to estimate the covariance matrix.
#### ctlruns2:   control runs used to estimate the covariance matrix.
#### ni:         number of locations for the observed responses.
#### C:          number of time steps for the observed responses.
#### nB:         number of bootstrap replications.
#### conf.level: confidence level.
#### cal_a:      indicator for calculating the corresponding a value in EE method.
eefp <- function(Xt, Y, m, ctlruns1, ctlruns2, ni, C,
                 nB = 0, conf.level = 0.9, cal_a = T) {
  p <- dim(Xt)[2]  ## number of forcings
  n <- length(Y)   ## number of forcings
  ## check number of time steps, if C > 1, use the pooled data to estimate
  ## the covariance matrix as a block diagonal matrix (a spatial-temporal covariance 
  ## matrix with temporal independence).
  if(C == 1){
    Sig <- lwRegcov(ctlruns1)
  } else {
    if(ni == 1) {
      Sig <- diag(C) %x% cov(t(poolDat(X = t(ctlruns1), C = C, ni = ni)))
    } else {
      Sig <- diag(C) %x% lwRegcov(t(poolDat(X = t(ctlruns1), C = C, ni = ni)))
    }
  }
  ## construct a diagonal matrix with each element being 1 / number of ensembles
  if(length(m) == 1){
    mv <- 1/m
  } else {
    mv <- diag(1/m)
  }
  X.o <- as.matrix(Xt)
  Y.o <- as.matrix(Y)
  ## separate the signal matrix
  X <- NULL
  for(i in 1:p){
    X <- cbind(X, matrix(X.o[, i], nrow = ni, ncol = C))
  }
  ## compute the value for estimating equation
  dotG_t <- lapply(1:C, function(i){
    t(t(X[, seq(i, (p - 1) * C + i, C), drop = FALSE]) %*% 
        solve(Sig[((i - 1) * ni + 1):(i * ni), ((i - 1) * ni + 1):(i * ni), drop = FALSE]) %*% 
        X[, seq(i, (p - 1)*C + i, C)] -
        ni * mv) %>% as.matrix})
  dotG <- Reduce(`+`, dotG_t)
  A <- solve(dotG)
  ## Point Estimation
  beta.hat <- A %*%
    Reduce(`+`, lapply(1:C, function(i){
      t(X[, seq(i, (p - 1) * C + i, C), drop = FALSE]) %*% 
        solve(Sig[((i - 1) * ni + 1):(i * ni), ((i - 1) * ni + 1):(i * ni), drop = FALSE]) %*% 
        Y[((i - 1)*ni + 1):(i*ni)] %>% as.matrix}))
  ## Residual and a
  Sig2 <- extract_block_diag(Sig, ni, C)
  res <- matrix(Y - X.o %*% beta.hat, nrow = ni, ncol = C)
  res_tmp2 <- sapply(1:C, function(i) inv_sr(Sig2[[i]]) %*% res[, i, drop = FALSE])
  res_tmp2 <- matrix(res_tmp2, ncol = C)  ## make sure it is a matrix
  res_tmp <- as.vector(res_tmp2)
  ## estimate alpha if required
  if(cal_a){
    if(length(m) > 1){
      a <- c(1/(var(res_tmp) - diag(mv) %*% beta.hat^2))
    } else {
      a <- c(1/(var(res_tmp) - mv * c(beta.hat)^2))
    }
  } else {
    a <- 1
  }
  ## pre-whitened residuals
  if(length(m) > 1) {
    res_pwh <- c(sqrt(a/(1 + a * diag(mv) %*% beta.hat^2))) * res_tmp
  } else {
    res_pwh <- sqrt(a/(1 + a * mv * c(beta.hat)^2)) * res_tmp
  }
  ## Variance and Confidence Interval
  l2 <- nrow(ctlruns2)
  ## define G function
  G_fun <- function(ep){
    G <- sapply(1:C, function(i){
      t(t(X[, seq(i, (p - 1)*C + i, C), drop = FALSE]) %*% 
          solve(Sig[((i - 1)*ni + 1):(i*ni), ((i - 1)*ni + 1):(i*ni)]) %*% 
          ep[((i - 1)*ni + 1):(i*ni)] + ni*mv%*%beta.hat) %>% as.matrix})
    ## solve the one demensional problem
    G <- matrix(G, nrow = p, ncol = C)
    return(G)
  }
  if(p > 1){
    Gsb <- sapply(1:l2, function(x){
      rowSums(G_fun(ctlruns2[x, ]))
    })
    Bsb <- c(1/a + diag(mv) %*% beta.hat^2) * cov(t(Gsb))
  } else {
    Gsb <- sapply(1:l2, function(x){
      sum(G_fun(ctlruns2[x, ]))
    })
    Bsb <- c(1/a + mv * beta.hat^2) * var(Gsb)
  }
  ## bootstrap to calculate the Variance if nB > 0
  Bb <- Bsb
  if(nB > 0) {
    l <- 2
    mblk <- sapply(1:nB,
                   function(x) {
                     idx1 <- sample(1: (C - l + 1), ceil(C / l), replace = TRUE)
                     idx2 <- c(outer(0:(l - 1), idx1, "+"))
                     if(length(idx2) > C){
                       idx2 <- idx2[-length(idx2)]
                     }
                     idx3 <- -c(outer((ni - 1):0, idx2*ni, "-"))
                     return(idx3)
                   })
    Gb <- apply(mblk, 2,
                function(x) {
                  rowSums(G_fun(Y.o[x,] - X.o[x, ] %*% beta.hat))
                })
    Gb <- matrix(Gb, nrow = p, ncol = nB)
    Bb <- cov(t(Gb))
  }
  norm.crt <- qnorm(1 - (1 - conf.level)/2)  ## critical value for normal approximation
  ## bootstrapped results
  Vsb <- A  %*% Bb %*% A
  sd_sb <- sqrt(diag(Vsb))
  # sd_sb <- sqrt(diag(A  %*% Bb %*% A))
  ci_sb <- cbind(beta.hat - norm.crt * sd_sb,
                 beta.hat + norm.crt * sd_sb)
  
  ## approximated variance and confidence interval estimations
  Vb <- A  %*% Bsb %*% A
  sd_b <- sqrt(diag(Vb))
  # sd_b <- sqrt(diag(A  %*% Bsb %*% A))
  ci_b <- cbind(beta.hat - norm.crt * sd_b,
                beta.hat + norm.crt * sd_b)
  
  ## calculate CI of a
  if(length(m) > 1){
    dotQ_t <- lapply(1:C, function(i){
      t(res_tmp2[, i]) %*% res_tmp2[, i] -
        (ni- 1/C)*t(beta.hat)%*% mv %*% beta.hat})
  } else {
    dotQ_t <- lapply(1:C, function(i){
      t(res_tmp2[, i])%*%res_tmp2[, i] -
        (ni- 1/C)*beta.hat^2*mv})
  }
  dotQ <- c(Reduce(`+`, dotQ_t))
  dC <- 1/dotQ
  Q_fun <- function(ep){
    Q <- sapply(1:C, function(i){
      a*t(ep[((i - 1) * ni + 1):(i * ni)]) %*% 
        solve(Sig[((i - 1) * ni + 1):(i * ni), ((i - 1) * ni + 1):(i * ni)]) %*% 
        ep[((i - 1) * ni + 1):(i * ni)]})
    return(Q)
  }
  Qsb <- sapply(1:l2, function(x){
    sum(Q_fun(ctlruns2[x, ])) 
  })
  if(length(m) > 1){
    Dsb <- c(1/a + diag(mv) %*% beta.hat^2) * var(Qsb)/a
    C2_t <- lapply(1:C, function(i){
      t(t(X[, seq(i, (p - 1)*C + i, C), drop = FALSE]) %*% 
          solve(Sig[((i - 1)*ni + 1):(i*ni), ((i - 1)*ni + 1):(i * ni)]) %*% 
          res[, i] + (ni- 1/C)* mv %*% beta.hat) %>% as.matrix})
  } else {
    Dsb <- c(1 / a + mv * beta.hat^2) * var(Qsb) / a
    C2_t <- lapply(1:C, function(i){
      t(t(X[, seq(i, (p - 1)*C + i, C), drop = FALSE]) %*% 
          solve(Sig[((i - 1)*ni + 1):(i*ni), ((i - 1) * ni + 1):(i * ni)]) %*% 
          res[, i] + (ni- 1/C) * mv * beta.hat) %>% as.matrix})
  }
  C2 <- 2 * a * c(Reduce(`+`, C2_t))
  ## variance of a
  sd_a <- sqrt(c(dC^2 * t(C2) %*% Vsb %*% C2) + dC^2 * Dsb)
  Z_stat <- c((a - 1) / sd_a)
  p_test_a <- 2 * min(pnorm(Z_stat), 1 - pnorm(Z_stat))
  ## confidence interval of a
  ci_a <- cbind(a - norm.crt * sd_a,
                a + norm.crt * sd_a)
  ## return the point estimate and confidence interval results as a list
  numbeta <- length(beta.hat)
  Xlab <- paste0("X", 1:numbeta)
  if (nB > 0){
    result <- list(beta = as.vector(beta.hat),  ## point estimate
                   var = Vb, ci = ci_b,  ## variance and CI
                   var.boot = Vsb, ci.boot = ci_sb,  ## bootstrapped variance and CI
                   a = a, sd.a = sd_a, ci.a = ci_a,  ## p-value and ci for testing alpha = 1
                   p.value.a = p_test_a,
                   residuals = c(res), residuals.prewhit = res_pwh)  ## residuals and pre-whitened residuals
    rownames(result$var.boot) <- colnames(result$var.boot) <- Xlab
    rownames(result$ci.boot) <- Xlab
    colnames(result$ci.boot) <- c("Lower bound", "Upper bound")
  } else {
    result <- list(beta = as.vector(beta.hat),
                   var = Vb, ci = ci_b,  ## approximated var and CI
                   a = a, sd.a = sd_a, ci.a = ci_a,
                   p.value.a = p_test_a,
                   residuals = c(res), residuals.prewhit = res_pwh)
  }
  ## rename the results
  names(result$beta) <- Xlab
  rownames(result$var) <- colnames(result$var) <- Xlab
  rownames(result$ci) <- Xlab
  colnames(result$ci) <- colnames(result$ci.a) <- c("Lower bound", "Upper bound")
  ## return the results
  return(result)
}

## fingerprinting model by EE method with missing value
#### Xt:         n \times p matrix, the ensemble model simulations, each column represents a
####             external forcing.
#### Y:          n \times 1 matrix, length ni (number of locations) x C (number of time steps), 
####             the corresponding observed values.
#### m:          p vector, the number of ensemble runs to obtain the signals of external forcings.
#### ctlruns1:   control runs used to estimate the covariance matrix.
#### ctlruns2:   control runs used to estimate the covariance matrix.
#### ni:         number of locations for the observed responses.
#### C:          number of time steps for the observed responses.
#### nB:         number of bootstrap replications.
#### conf.level: confidence level.
#### ridge:      coefficient for the shrinkage to handle missing value.

## check the missing pattern
eefp_mis <- function(Xt, Y, m, ctlruns1, ctlruns2, ni, C, 
                     ridge = 0, nB = 0, conf.level = 0.9) {
  ## check the missing value
  mis <- which(is.na(Y))
  ## mask the missing value in the X and covariance matrix
  Xt[mis, ] <- NA
  p <- dim(Xt)[2]
  n <- length(Y)
  l2 <- dim(ctlruns2)[1]
  
  if(C == 1){
    Sig <- lwRegcov(ctlruns1)
  } else {
    if(ni == 1) {
      Sig <- diag(C) %x% cov(t(poolDat(X = t(ctlruns1), C = C, ni = ni)))
    } else {
      Sig <- diag(C) %x% lwRegcov(t(poolDat(X = t(ctlruns1), C = C, ni = ni)))
    }
  }
  Sig[mis, ] <- NA
  Sig[, mis] <- NA
  if(length(m) == 1){
    mv <- 1/m
  } else {
    mv <- diag(1/m)
  }
  X.o <- as.matrix(Xt)
  Y.o <- as.matrix(Y)
  X <- NULL
  for(i in 1:p){
    X <- cbind(X, matrix(X.o[, i], nrow = ni, ncol = C))
  }
  ## compute the value for estimating equation
  dotG_t <- lapply(1:C, function(i){
    t(t(remove_empty(which = c("rows", "cols"), dat = X[, seq(i, (p - 1) * C + i, C), drop = FALSE])) %*% 
        solve(remove_empty(which = c("rows", "cols"), dat = Sig[((i - 1) * ni + 1):(i*ni), ((i - 1) * ni + 1):(i * ni), drop = FALSE])) %*%
        remove_empty(which = c("rows", "cols"), dat = X[, seq(i, (p - 1) * C + i, C), drop = FALSE]) -
        ni * mv) %>% as.matrix})
  dotG <- Reduce(`+`, dotG_t)
  A <- solve(dotG + ridge * diag(p))
  beta.hat <- A %*%
    Reduce(`+`, lapply(1:C, function(i){
      t(remove_empty(which = c("rows", "cols"), dat = X[, seq(i, (p - 1) * C + i, C), drop = FALSE])) %*% 
        solve(remove_empty(which = c("rows", "cols"), dat = Sig[((i - 1) * ni + 1):(i * ni), ((i - 1) * ni + 1):(i * ni), drop = FALSE])) %*% 
        remove_empty(which = c("rows", "cols"), dat = matrix(Y[((i - 1) * ni + 1):(i * ni)])) %>% as.matrix}))
  G_fun <- function(ep){
    ep[mis] <- NA
    G <- sapply(1:C, function(i){
      t(t(remove_empty(which = c("rows", "cols"), dat = X[, seq(i, (p - 1) * C + i, C), drop = FALSE])) %*% 
          solve(remove_empty(which = c("rows", "cols"), dat = Sig[((i - 1) * ni + 1):(i * ni), ((i - 1) * ni + 1):(i * ni), drop = FALSE])) %*% 
          remove_empty(which = c("rows", "cols"), dat = matrix(ep[((i - 1) * ni + 1):(i * ni)])) +
          ni * mv %*% beta.hat) %>% as.matrix})
    ## solve the one demensional problem
    G <- matrix(G, nrow = p, ncol = C)
    return(G)
  }
  Gsb <- sapply(1:l2, function(x){
    rowSums(G_fun(ctlruns2[x, ]))
  })
  Bsb <- c(1 + diag(mv) %*% beta.hat^2) * cov(t(Gsb))
  ## bootstrap to calculate the Variance if nB > 0
  Bb <- Bsb
  if(nB > 0) {
    l <- 2
    mblk <- sapply(1:nB,
                   function(x) {
                     idx1 <- sample(1: (C - l + 1), ceil(C/l), replace = TRUE)
                     idx2 <- c(outer(0:(l - 1), idx1, "+"))
                     if(length(idx2) > C){
                       idx2 <- idx2[-length(idx2)]
                     }
                     idx3 <- -c(outer((ni - 1):0, idx2*ni, "-"))
                     return(idx3)
                   })
    Gb <- apply(mblk, 2,
                function(x) {
                  rowSums(G_fun(Y.o[x, ] - X.o[x, ]%*%beta.hat))
                })
    ## fixed nB issue
    Gb <- matrix(Gb, nrow = p, ncol = nB)
    Bb <- cov(t(Gb))
  }
  norm.crt <- qnorm(1 - (1 - conf.level)/2)  ## critical value for normal approximation
  ## bootstrapped results
  Vsb <- A  %*% Bb %*% A
  sd_sb <- sqrt(diag(Vsb))
  # sd_sb <- sqrt(diag(A  %*% Bb %*% A))
  ci_sb <- cbind(beta.hat - norm.crt * sd_sb,
                 beta.hat + norm.crt * sd_sb)
  
  ## approximated variance and confidence interval estimations
  Vb <- A  %*% Bsb %*% A
  sd_b <- sqrt(diag(Vb))
  # sd_b <- sqrt(diag(A  %*% Bsb %*% A))
  ci_b <- cbind(beta.hat - norm.crt * sd_b,
                beta.hat + norm.crt * sd_b)
  ## summarize the results as a list
  numbeta <- length(beta.hat)
  Xlab <- paste0("X", 1:numbeta)
  if (nB > 0){
    result <- list(beta = as.vector(beta.hat),  ## point estimate
                   var = Vb, ci = ci_b,  ## variance and CI
                   var.boot = Vsb, ci.boot = ci_sb)  ## bootstrapped variance and CI
    rownames(result$var.boot) <- colnames(result$var.boot) <- Xlab
    rownames(result$ci.boot) <- Xlab
    colnames(result$ci.boot) <- c("Lower bound", "Upper bound")
  } else {
    result <- list(beta = as.vector(beta.hat),
                   var = Vb, ci = ci_b)
  }
  ## rename the results
  names(result$beta) <- Xlab
  rownames(result$var) <- colnames(result$var) <- Xlab
  rownames(result$ci) <- Xlab
  colnames(result$ci) <- c("Lower bound", "Upper bound")
  return(result)
}
## internal utility functions

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
  sample.cov <- cov(X)
  Ip <- diag(p)
  m <- sum(diag(sample.cov %*% Ip)) / p ## first estimate in L&W
  Xp <- sample.cov - m * Ip
  d2 <- sum(diag(Xp %*% t(Xp))) / p  ## second estimate in L&W
  
  bt <- (diag(X %*% t(X))^2 - 2 * diag(X %*% sample.cov %*% t(X)) + rep(1, n) *
           sum(diag(sample.cov %*% sample.cov))) / p
  
  bb2 <- 1 / n^2 * sum(bt)
  b2 <- min(bb2,d2)  ## third estimate in L&W
  a2 <- d2 - b2  ## fourth estimate in L&W
  ## the regularized estimate of the covariance matrix
  ## b2 * m / d2 * Ip + a2 / d2 * sample.cov
  eigen(sample.cov)$vector %*% 
    diag(b2 * m / d2 + a2 / d2 * eigen(sample.cov)$values) %*% 
    t(eigen(sample.cov)$vector)
}

# pool data
poolDat <- function(X, C, ni){
  X.pool <- X[1:ni, , drop = FALSE]
  for(i in 2:C){
    X.pool <- cbind(X.pool, 
                    X[((i - 1)*ni + 1):(i*ni),  , drop = FALSE])
  }
  return(X.pool)
}

## for the residual and a values
# extract block diagonal matrix
extract_block_diag <- function(x, ni, C){
  blk <- list()
  for(i in 1:C){
    blk_tmp <- x[((i - 1)*ni + 1):(ni*i), ((i - 1)*ni + 1):(ni*i)]
    if(i == 1){
      blk <- list(as.matrix(blk_tmp))
    } else{
      # blk <- list(blk, blk_tmp)
      blk <- append(blk, list(as.matrix(blk_tmp)))
    }
  }
  return(blk)
}

# calculate inverse matrix
inv_sr <- function(x){
  if(length(x)  == 1) {
    inv_sr <- 1 / sqrt(x)
  } else {
    ei <- eigen(x)
    eive <- ei$vectors
    inv_sr <- eive %*% diag(1 / sqrt(ei$values)) %*% t(eive)
  }
  return(inv_sr)
}

# prewhitening
pre_wh_ctlr <- function(vec1, Sig, ni, C){
  vec2 <- matrix(vec1, nrow = ni, ncol = C)
  Sig2 <- extract_block_diag(Sig, ni, C)
  rt <- sapply(1:C, function(i) inv_sr(Sig2[[i]])%*%vec2[, i])
  # rt <- inv_sr(Sig) %*% vec1
  return(as.vector(rt))
}


#### create a block Toeplitz covariance matrix from control runs
lwRegcov5 <- function(X, C, ni) {
  sample.cov <- cov(X)
  ## apply averages on correlation matrix to form a toeplitz block might help.  
  sample.cov <- toe_cov(sample.cov, C, ni)
  n <- nrow(X)
  p <- ncol(X)
  Ip <- diag(p)
  m <- sum(diag(sample.cov %*% Ip)) / p ## first estimate in L&W
  Xp <- sample.cov - m * Ip
  d2 <- sum(diag(Xp %*% t(Xp))) / p  ## second estimate in L&W
  
  bt <- (diag(X %*% t(X))^2 - 2 * diag(X %*% sample.cov %*% t(X)) + rep(1, n) *
           sum(diag(sample.cov %*% sample.cov))) / p
  
  bb2 <- 1 / n^2 * sum(bt)
  b2 <- min(bb2,d2)  ## third estimate in L&W
  a2 <- d2 - b2  ## fourth estimate in L&W
  ## the regularized estimate of the covariance matrix
  eg <- eigen(sample.cov) ## eigen is time consuming; only need it once
  lambda <- eg$value
  thres <- 1e-5
  lambda[lambda < thres] <- 0
  Shalf <- sqrt(b2 * m / d2 + a2 / d2 * lambda) * t(eg$vector)
  S <- crossprod(Shalf)
  return(S)
}

#### compute the covariance matrix with pooled data
ave_cov <- function(S, C, ni){
  S.pool <- poolDat(t(S), C, ni)
  rt <- cov(t(S.pool))
}

#### Toeplitz 
#### a function to make a symmetric matrix to be a block Toeplitz matrix
toe_cov <- function(S, C, ni){
  Si.tmp <- replicate(C, matrix(0, nrow = ni, ncol = ni), simplify = FALSE)
  S.tmp <- S
  for(i in 1:C){
    Si.tmp[[1]] <- Si.tmp[[1]] + S.tmp[((i - 1)*ni + 1):(i*ni), ((i - 1)*ni + 1):(i*ni)]
  }
  for(j in 2:C) {
    S.tmp <- S.tmp[-(1:ni), -((dim(S.tmp)[1] - ni +1):dim(S.tmp)[1])]
    for(i in 1:(C - j + 1)){
      Si.tmp[[j]] <- Si.tmp[[j]] + S.tmp[((i - 1)*ni + 1):(i*ni), ((i - 1)*ni + 1):(i*ni)]
    }
  }
  for(k in 1:C){
    Si.tmp[[k]] <- Si.tmp[[k]]/(C - k + 1)
  }
  # impose decreasing property
  for(k in 3:C){
    idx <- Si.tmp[[k]] > Si.tmp[[k - 1]]
    if(sum(idx) > 0){
      Si.tmp[[k]][idx] <- Si.tmp[[k - 1]][idx]
    }
  }
  rt <- matrix(data = 0, nrow = C*ni, ncol = C*ni)
  Si.tmp.t <- Map("t", Si.tmp)
  rt[1:ni, ] <- do.call(cbind, Si.tmp)
  rt[, 1:ni] <- do.call(rbind, Si.tmp.t)
  for(l in 2:C) {
    rt[((l - 1)*ni + 1):(l*ni), (ni + 1):(C*ni)] <- rt[((l - 2)*ni + 1):((l - 1)*ni), 1:((C - 1)*ni)]
  }
  return(rt)
}
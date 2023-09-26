##' Regularized estimators for covariance  matrix.
##'
##' This function estimate the covariance matrix under l2 loss and 
##' minimum variance loss, provide linear shrinkage estimator under 
##' l2 loss and nonlinear shrinkage estimator under minimum variance loss.
##' 
##' @param Z n*p matirx with sample size n and dimension p. 
##' Replicates for computing the covariance matrix, should be centered. 
##' @param method methods used for estimating the covariance matrix.
##' @param bandwidth bandwidth for the "mv" estimator, 
##' default value are set to be list in (0.2, 0.5).
##' @return regularized estimate of covariance matrix.
##' @author Yan Li
##' @keywords htest
##' @references \itemize{ \item  Olivier Ledoit and Michael Wolf,
##' A well-conditioned estimator for large-dimensional 
##' covariance matrices, 2004, JMA. 
##' \item Olivier Ledoit and Michael Wolf, 
##' Direct nonlinear shrinkage estimation of 
##' large-dimensional covariance matrices, Working Paper No. 264, UZH.
##' \item Li et al,
##' Regularized fingerprinting in detection and attribution of climate change 
##' with weight matrix optimizing the efficiency in 
##' scaling factor estimation, 2023, The Annals of Applied Statistics.}
##' @examples
##' ## randomly generate a n * p matrix where n = 50, p = 100
##' Z <- matrix(rnorm(50 * 100), nrow = 50, 100)
##' ## linear shrinkage estimator under l2 loss
##' Cov.est <- Covest(Z, method = "l2")$output
##' ## nonlinear shrinkage estimator under minimum variance loss
##' Cov.est <- Covest(Z, method = "mv", bandwidth = 0.35)$output
##' @importFrom Iso pava
##' @export Covest
Covest <- function(Z, method = c("mv", "l2"), bandwidth = NULL) {
  method <- match.arg(method)
  ## create initial path for the bandwidth
  if (is.null(bandwidth)) {
    bandwidth <- seq(0.2, 0.5, by = 0.025)
  }
  
  if (!method %in% c("l2", "stein", "mv")) {
    stop("please select one objective function from 'l2', 'stein' or 'mv'")
  }
  output <- NULL
  if (method == "l2") {
    ## Linear shrinkage estimator under frobenious loss function
    output <- lwRegcov(Z)
    h.par <- 0
  } else if (method == "mv") {
    if(length(bandwidth) > 1) {
      Z.cv <- cvFolds(Z[1:(floor(nrow(Z) / 5) * 5), ])
      SD <- lapply(Z.cv$test, cov)
      ## Nonlinear shrinkage estimator under minimum variance loss function
      mv.select <-NULL
      for (h.par in bandwidth) {
        output.temp <- lapply(Z.cv$train, 
                              function(x) {
                                lwMvcov(x, h.par)
                              })
        mv.select <- c(mv.select, mean(sapply(1:length(SD), 
                                              function(ind) {
                                                mvDist(output.temp[[ind]], 
                                                       SD[[ind]])
                                              })))
      }
      h.par <- bandwidth[which.min(mv.select)]
      output <- lwMvcov(Z, h.par)
    } else {
      h.par <- bandwidth
      output <- lwMvcov(Z, h.par)
    }
  } else {
    ## Nonlinear shrinkage estimator under stein's loss function
    ## output <- lwSteincov(Z)
    stop("Unknown type of estimator")
  }
  list(output = output, bandwith = h.par)
}

mvDist <- function(Cov.hat, Cov) {
  Cov.hati <- solve(Cov.hat)
  Min <- (sum(diag(Cov.hati %*% Cov %*% Cov.hati)) / nrow(Cov)) / 
    (sum(diag(Cov.hati)) / nrow(Cov))^2
  Min
}

cvFolds <- function(input, k = 5) {
  ## input should be a n*p matrix with n observations and p covariates
  ## number of observations
  n <- nrow(input)
  n.set <- floor(n / k)
  if (n.set < n /k) {
    stop("number of observations in each group must be equal")
  }
  lower <- seq(1, n, by = n.set)
  upper <- seq(n.set, n, by = n.set)
  train <- test <- NULL
  for(i in 1:k) {
    test <- c(test, list(input[c(lower[i]:upper[i]), ]))
    train <- c(train, list(input[-c(lower[i]:upper[i]), ]))
  }
  list(test=test, train=train)
}


## compute the regularized estimate of the covariance matrix
## References:
## Ledoit O., M. Wolf (2004) A well-conditioned estimator for
## large-dimensional covariance matrices.
## Journal of Multivariate Analysis, 88(2): 365-411.
## L&W estimate
lwRegcov <- function(Z) {
  ## input:
  ##   the independent replicates of the random vector, denoted by matrix Z
  ##   each column represents a random variable. Thus Z is n*p matrix
  ## output:
  ##   regularized covariance matrix
  n <- nrow(Z)
  p <- ncol(Z)
  sample.cov <- cov(Z)
  Ip <- diag(p)
  m <- sum(diag(sample.cov %*% Ip)) / p ## first estimate in L&W
  Zp <- sample.cov - m * Ip
  d2 <- sum(diag(Zp %*% t(Zp))) / p  ## second estimate in L&W
  
  bt <- (diag(Z %*% t(Z))^2 - 2 * diag(Z %*% sample.cov %*% t(Z)) + rep(1, n) *
           sum(diag(sample.cov %*% sample.cov))) / p
  
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
lwMvcov <- function(Z, h.par) {
  ## input:
  ##   the sample replicates of the random vector,
  ##   denoted by n*p matrix Z.
  ## output:
  ##   nonlinear shrinkage covariance matrix under minimum variance
  ##   loss function.
  if(max(abs(colMeans(Z))) < 1e-12) {
    n <- nrow(Z) - 1
  } else {
    n <- nrow(Z)
  }
  p <- ncol(Z)
  ## sample covariance matrix
  sample.cov <- cov(Z)
  ## sample.cov <- t(Z) %*% Z / n
  ## eigen value and vectors
  lambda <- eigen(sample.cov)$value
  eigen.v <- eigen(sample.cov)$vectors
  temOd <- order(lambda)
  lambda <- lambda[temOd]
  eigen.v <- eigen.v[, temOd]
  ## check whether the eigenvalues are close to zero
  ## avoid effects of extremely small eigenvalues
  test <- sum(sapply(lambda, 
                     function(x) {
                       x < 5 * 10^{-5}
                     }) != "TRUE")
  if(test < p) {
    n <- test
  }
  if(p == n) {
    n <- n - 1
  }
  ## compute the direct kernal estimator
  lambda.tilde <- lambda[max(1, p - n  + 1):p]
  h <- n^(-h.par)
  L <- matrix(rep(lambda.tilde, min(p, n)), length(lambda.tilde), min(p,n))
  TemM <- 4 * t(L)^2 * h^2 - (L - t(L))^2
  TemM[TemM < 0] = 0
  ftilde <- rowMeans(sqrt(TemM) / (2 * pi * h^2 * t(L)^2))
  TemM <- -(4 * t(L)^2 * h^2 - (L - t(L))^2)
  TemM[TemM < 0] = 0
  Hftilde <- rowMeans((sign(L - t(L)) * sqrt(TemM) - L + t(L)) /
                        (2 * pi * h^2 * t(L)^2))
  if(p <= n) {
    ## non-singular case
    dtilde <- lambda.tilde / ((pi * p / n * lambda.tilde * ftilde)^2 +
                                (1 - p/n - pi * p/n * lambda.tilde * Hftilde)^2)
  } else {
    ## singular case
    Hftilde0 <- (1 - sqrt(1 - 4 * h^2)) / (2 * pi * h^2) * mean(1 / lambda.tilde)
    dtilde0 <- 1 / (pi * (p - n) / n * Hftilde0)
    dtilde1 <- lambda.tilde / (pi^2 * lambda.tilde^2 * (ftilde^2 + Hftilde^2))
    dtilde <- c(rep(dtilde0, p - n), dtilde1)
  }
  ## Pool Adjacent Violators algorithm.
  dtilde <- pava(dtilde)
  ## output the covariance matrix
  eigen.v %*% diag(dtilde) %*% t(eigen.v)
}

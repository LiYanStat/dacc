##' Regularized optimal fingerprinting with linearly optimal weight matrix.
##'
##' Implements the regularized fingerprinting method proposed by Li, Chen, and
##' Yan (2025). The method sweeps over a grid of regularization parameters
##' \eqn{\lambda} for the weight matrix \eqn{(\hat{\Sigma} + \lambda I)^{-1}}
##' and selects the value that minimizes the total asymptotic estimation
##' uncertainty across all scaling factors. Variance estimation uses random
##' matrix theory corrections that account for the high-dimensional nature of
##' climate covariance estimation.
##'
##' @param X an \eqn{N \times p} matrix of model-simulated responses to
##'   \eqn{p} external forcings (fingerprint patterns). Each column corresponds
##'   to one forcing.
##' @param Y an \eqn{N \times 1} matrix of observed climate variables.
##' @param nruns.X a numeric vector of length \eqn{p} giving the number of
##'   ensemble runs used to estimate each column of \code{X}.
##' @param ctlruns an \eqn{m \times N} matrix of pre-industrial control run
##'   segments used to estimate the covariance structure of internal climate
##'   variability.
##' @param conf.level confidence level for the scaling factor intervals.
##'   Default is \code{0.90}.
##' @param lambda an optional numeric vector of regularization parameters to
##'   evaluate. If \code{NULL} (default), an automatic grid of 20 values
##'   spanning \eqn{[0.01\bar{\tau},\; 10\bar{\tau}]} is constructed, where
##'   \eqn{\bar{\tau} = N^{-1}\,\mathrm{tr}(\hat{S})} is the mean eigenvalue of
##'   the sample covariance of the control runs (arXiv:2505.04070, §2.3).
##'
##' @return A list with the following components:
##'   \itemize{
##'     \item \code{beta.list}: \eqn{p \times m} matrix of scaling factor
##'       estimates for each value of \eqn{\lambda}.
##'     \item \code{var.list}: \eqn{p \times p \times m} array of
##'       variance-covariance matrices for each \eqn{\lambda}.
##'     \item \code{ci.list}: \eqn{2 \times p \times m} array of confidence
##'       intervals for each \eqn{\lambda}. The first dimension indexes the
##'       lower and upper bounds.
##'     \item \code{sd.list}: \eqn{p \times m} matrix of standard deviations
##'       for each \eqn{\lambda}.
##'     \item \code{beta.optim}: length-\eqn{p} vector of scaling factor
##'       estimates at the optimal \eqn{\lambda}.
##'     \item \code{var.optim}: \eqn{p \times p} variance-covariance matrix
##'       at the optimal \eqn{\lambda}.
##'     \item \code{ci.optim}: \eqn{2 \times p} matrix of confidence interval
##'       bounds at the optimal \eqn{\lambda}.
##'     \item \code{sd.optim}: length-\eqn{p} vector of standard deviations
##'       at the optimal \eqn{\lambda}.
##'     \item \code{lambda}: numeric vector of all \eqn{\lambda} values
##'       evaluated.
##'     \item \code{Xi.optim}: \eqn{p \times p} matrix \eqn{N \hat{\Xi}(\hat{\lambda}_{\rm opt})},
##'       the \eqn{N}-scaled asymptotic covariance as defined in Theorem 1 of
##'       arXiv:2505.04070. Equals \code{var.optim * N}.
##'     \item \code{joint.cr}: a list giving the joint confidence region for
##'       \eqn{\beta} (arXiv:2505.04070, §2.3). Contains:
##'       \code{stat} (the test statistic
##'       \eqn{\hat{\beta}^{\top}\hat{\Xi}^{-1}(\hat{\lambda}_{\rm opt})\hat{\beta}}),
##'       \code{bound} (the \eqn{\chi^2_p(1-\alpha)} critical value), and
##'       \code{detected} (\code{TRUE} if the test statistic exceeds the bound,
##'       i.e., the zero vector lies outside the joint confidence region).
##'   }
##'
##' @author Yan Li, Kun Chen, Jun Yan
##' @keywords climate detection attribution fingerprinting regularization
##'
##' @references
##' \itemize{
##'   \item Li, Y., Chen, K., Yan, J. (2025). Regularized Fingerprinting with
##'     Linearly Optimal Weight Matrix in Detection and Attribution of Climate
##'     Change. \emph{Journal of Climate}.
##'   \item Ledoit, O., Wolf, M. (2004). A well-conditioned estimator for
##'     large-dimensional covariance matrices. \emph{Journal of Multivariate
##'     Analysis}, 88(2), 365--411.
##' }
##'
##' @examples
##' data(globalDat)
##' nomis <- which(!is.na(globalDat$Y))                                                                                                                                 
##' result <- fingerprintOPT(                                                                                                                                         
##'   X       = globalDat$Xtilde[nomis, ],
##'   Y       = globalDat$Y[nomis, , drop = FALSE],                                                                                                                     
##'   nruns.X = globalDat$nruns.X[c("ANT", "NAT")],
##'   ctlruns = globalDat$ctlruns[, nomis])
##' result$beta.optim
##' result$ci.optim
##'
##' @export fingerprintOPT
##' @importFrom stats cov qchisq qnorm
fingerprintOPT <- function(X, Y, nruns.X, ctlruns,
                           conf.level = 0.90,
                           lambda = NULL) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("forcing ", 1:ncol(X))
  }
  ## set the range of lambda
  if (is.null(lambda)) {
    ## arXiv:2505.04070, §2.3: grid anchored to tau-bar = mean eigenvalue of S
    S.hat   <- cov(ctlruns)
    tau.bar <- mean(eigen(S.hat, only.values = TRUE)$values)
    lambda  <- seq(0.01 * tau.bar, 10 * tau.bar, length = 20)
  }
  ## fit the TLS model across the full lambda grid
  output <- tlsLm.opt(X = X, Y = Y,
                      nruns.X = nruns.X,
                      conf.level = conf.level,
                      lambda = lambda,
                      ctlruns = ctlruns)
  ## select the optimal lambda: minimises total squared standard deviation
  ind.opt <- which.min(colSums(output$sd^2))
  output$beta.optim <- output$beta.list[, ind.opt]
  output$var.optim  <- output$var.list[, , ind.opt]
  output$ci.optim   <- output$ci.list[, , ind.opt]
  output$sd.optim   <- output$sd[, ind.opt]
  output$lambda     <- lambda
  ## arXiv:2505.04070, Theorem 1: N-scaled asymptotic covariance matrix
  output$Xi.optim <- output$var.optim * nrow(Y)
  ## arXiv:2505.04070, §2.3: joint confidence region for beta
  p <- length(output$beta.optim)
  jcr.stat <- as.numeric(t(output$beta.optim) %*%
                solve(output$var.optim) %*% output$beta.optim)
  jcr.bound <- qchisq(conf.level, df = p)
  output$joint.cr <- list(
    stat     = jcr.stat,
    bound    = jcr.bound,
    detected = jcr.stat >= jcr.bound
  )
  output
}


## Total least squares estimation across a grid of regularization parameters.
##
## Internal worker called by fingerprintOPT(). For each value of lambda,
## computes TLS scaling factor estimates and their asymptotic variances using
## RMT-corrected auxiliary quantities from estimators().
##
## @param X N x p matrix of (unnormalized) response patterns.
## @param Y N x 1 matrix of observations.
## @param nruns.X p-vector of ensemble sizes per forcing.
## @param conf.level confidence level for intervals.
## @param lambda numeric vector of regularization parameters.
## @param ctlruns m x N control run matrix.
## @return list with beta.list, var.list, ci.list, sd.
tlsLm.opt <- function(X, Y, nruns.X,
                      conf.level,
                      lambda,
                      ctlruns) {
  if (!is.matrix(Y)) {
    stop("Y should be an N x 1 matrix")
  }
  if (nrow(X) != nrow(Y)) {
    stop("X and Y must have the same number of rows")
  }
  n <- nrow(Y)         ## number of observation dimensions (S x T)
  p <- ncol(X)         ## number of forcings
  m <- length(lambda)  ## number of lambda values in the grid
  ## normalise X columns by sqrt(nruns.X) to equalise ensemble sampling
  ## variability across forcings before computing TLS
  Xt <- X * t(sqrt(nruns.X) %*% matrix(1, 1, n))
  if (length(nruns.X) != 1) {
    Dn.X <- diag(sqrt(nruns.X))
  } else {
    Dn.X <- sqrt(nruns.X)
  }
  ## compute RMT-based auxiliary variance arrays across all lambda values
  est.array <- estimators(Z = ctlruns,
                          lambda.set = lambda,
                          X.tilde = Xt,
                          n.i.set = rep(1, p))
  ## eigendecomposition of sample covariance from control runs
  sample.cov <- cov(ctlruns)
  tmp.eig    <- eigen(sample.cov)
  eig.val    <- tmp.eig$values
  eig.vec    <- tmp.eig$vectors
  ## TLS solver for a single lambda value:
  ## beta = (Xt'D Xt - lambda_min * I)^{-1} Xt'D Y
  ## where D = diag(1 / (eigenvalues + lambda)) and
  ## lambda_min is the smallest eigenvalue of the regularised augmented matrix
  tls.fit <- function(lam, Xt, Y, eig.val, eig.vec) {
    M      <- t(eig.vec) %*% cbind(Xt, Y)
    MtM    <- t(M) %*% diag(1 / (eig.val + lam)) %*% M
    lam.min <- tail(eigen(MtM)$values, 1)
    as.vector(solve(MtM[1:p, 1:p] - lam.min * diag(p)) %*% MtM[1:p, p + 1])
  }
  ## scaling factor estimates for every lambda in the grid (p x m matrix)
  beta.list1 <- sapply(lambda, tls.fit,
                       Xt = Xt, Y = Y,
                       eig.val = eig.val, eig.vec = eig.vec)
  ## asymptotic variance matrix for every lambda (p x p x m array)
  var.list1 <- sapply(seq_len(m), function(ilam) {
    (1 + sum(beta.list1[, ilam]^2)) *
      solve(est.array$Delta1[, , ilam]) %*%
      (est.array$Delta2[, , ilam] +
         est.array$K[ilam] *
         solve(diag(p) + beta.list1[, ilam] %*% t(beta.list1[, ilam]))) %*%
      solve(est.array$Delta1[, , ilam]) / n
  }, simplify = "array")
  ## back-transform by Dn.X to return estimates on the original scale
  beta.list <- Dn.X %*% beta.list1
  var.list  <- array(apply(var.list1, 3, function(v) Dn.X %*% v %*% Dn.X),
                     dim(var.list1))
  colnames(beta.list)     <- lambda
  dimnames(var.list)[[3]] <- lambda
  ## confidence intervals from asymptotic normality
  Z.crt  <- qnorm((1 - conf.level) / 2, lower.tail = FALSE)
  sd.mat <- sqrt(apply(var.list, 3, diag))
  ## Build ci as 2 x p x m: dim 1 = bound (1=lower, 2=upper),
  ## dim 2 = forcing, dim 3 = lambda.  Direct assignment avoids the
  ## column-major reshaping ambiguity of array(rbind(...)).
  ci <- array(NA_real_, dim = c(2, p, m))
  ci[1, , ] <- beta.list - Z.crt * sd.mat
  ci[2, , ] <- beta.list + Z.crt * sd.mat
  ## attach forcing names to all outputs
  forcing.names <- colnames(X)
  dimnames(beta.list)[[1]] <-
    dimnames(sd.mat)[[1]]  <-
    dimnames(var.list)[[1]] <-
    dimnames(var.list)[[2]] <-
    dimnames(ci)[[2]] <- forcing.names
  dimnames(ci)[[1]] <- c("lower", "upper")
  dimnames(ci)[[3]] <- as.character(lambda)
  list(beta.list = beta.list,
       var.list  = var.list,
       ci.list   = ci,
       sd        = sd.mat)
}


## Compute auxiliary matrices Delta1, Delta2, and K for asymptotic variance.
##
## Uses random matrix theory in the Marchenko-Pastur regime to correct for
## finite-sample bias introduced by high-dimensional covariance estimation.
## Called internally by tlsLm.opt().
##
## @param Z n.clt x N.dim control run matrix.
## @param lambda.set numeric vector of regularization parameters.
## @param X.tilde N.dim x p normalised response pattern matrix.
## @param n.i.set p-vector: reciprocal sampling ratios (all 1 in standard use).
## @return list with Delta1 (p x p x m), Delta2 (p x p x m), K (m-vector).
estimators <- function(Z, lambda.set, X.tilde, n.i.set) {
  n.clt <- nrow(Z)         ## number of control run segments
  N.dim <- ncol(Z)         ## observation dimension
  p.dim <- ncol(X.tilde)   ## number of forcings
  gamma <- N.dim / n.clt   ## aspect ratio
  m     <- length(lambda.set)
  D.sq  <- diag(1 / n.i.set)
  ## sample covariance and eigendecomposition
  S.hat   <- cov(Z)
  tmp.eig <- eigen(S.hat)
  eig.val <- tmp.eig$values
  U       <- tmp.eig$vectors
  ## Stieltjes transform quantities m1 and m2 for each lambda
  mF      <- sapply(lambda.set, m1, emp.eig = eig.val)
  mFprime <- sapply(lambda.set, m2, emp.eig = eig.val)
  ## Theta1 and Theta2: one scalar per lambda value
  Theta1 <- (1 - lambda.set * mF) /
    (1 - gamma * (1 - lambda.set * mF))
  Theta2 <- (1 + gamma * Theta1)^2 *
    (Theta1 - lambda.set * (mF - lambda.set * mFprime) /
       (1 - gamma * (1 - lambda.set * mF))^2)
  ## broadcast Theta1 and Theta2 into p x p x m arrays
  Theta1.arr <- sapply(Theta1, function(x) x * D.sq, simplify = "array")
  Theta2.arr <- sapply(Theta2, function(x) x * D.sq, simplify = "array")
  ## G1 and G2: p x p matrices for each lambda
  Q <- (1 / sqrt(N.dim)) * t(X.tilde) %*% U
  G1.arr <- sapply(lambda.set, function(lam) {
    Q %*% diag(1 / (eig.val + lam)) %*% t(Q)
  }, simplify = "array")
  G2.arr <- sapply(lambda.set, function(lam) {
    Q %*% diag(1 / (eig.val + lam)^2) %*% t(Q)
  }, simplify = "array")
  ## Delta1 = G1 - Theta1 * D_sq
  Delta1 <- G1.arr - Theta1.arr
  ## Delta2: enforce positive definiteness; fall back to tmp if violated
  tmp <- sweep(G1.arr, 3, (1 + gamma * Theta1)^2, "*") -
    sweep(G2.arr, 3, lambda.set * (1 + gamma * Theta1)^2, "*")
  Delta2 <- tmp - Theta2.arr
  for (i in seq_len(m)) {
    if (!all(eigen(Delta2[, , i])$values > 0)) {
      Delta2[, , i] <- tmp[, , i]
    }
  }
  list(Delta1 = Delta1, Delta2 = Delta2, K = Theta2)
}


## Stieltjes transform helper: mean(1 / (eigenvalues + lambda)).
m1 <- function(lambda, emp.eig) {
  mean(1 / (emp.eig + lambda))
}


## Stieltjes transform helper: mean(1 / (eigenvalues + lambda)^2).
m2 <- function(lambda, emp.eig) {
  mean(1 / (emp.eig + lambda)^2)
}

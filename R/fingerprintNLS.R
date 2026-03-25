##' Adaptable fingerprinting with nonlinear shrinkage.
##'
##' Implements the adaptable fingerprinting method proposed by Li and Li (2025).
##' The weight matrix is a cubic polynomial function of the sample covariance,
##' \eqn{W(f) = f(\hat{S})} where
##' \eqn{f(x) = 1 / (\ell_3 x^3 + \ell_2 x^2 + \ell_1 x + 1)},
##' applied eigenvalue-by-eigenvalue.  The three coefficients
##' \eqn{(\ell_3, \ell_2, \ell_1)} are selected by minimizing
##' \eqn{\mathrm{tr}(\hat{\Xi}(f))}, the total asymptotic estimation
##' uncertainty, over a grid of monotonically admissible candidates.
##' Variance estimation uses random matrix theory corrections that remain valid
##' in the high-dimensional regime \eqn{N / m \to \gamma \in (0, \infty)}.
##'
##' Two inflation factors \eqn{\sigma_X^2} and \eqn{\sigma_Z^2} are jointly
##' estimated from the raw ensemble members.  These account for the common
##' situation where climate model variability is proportionally mis-scaled
##' relative to the observations.  Setting \code{sigma.X2 = sigma.Z2 = 1}
##' recovers the special case of no inflation (the assumption made by
##' \code{fingerprintOPT}).
##'
##' @param X.tilde.list a list of \eqn{p} matrices, one per forcing.  Element
##'   \eqn{i} is an \eqn{n_i \times N} matrix whose rows are the individual
##'   ensemble members for forcing \eqn{i}.  The ensemble mean of each element
##'   is used as the fingerprint pattern.  The raw members are required for
##'   inflation factor estimation.
##' @param Y an \eqn{N \times 1} matrix of observed climate variables.
##' @param nruns.X a numeric vector of length \eqn{p} giving the number of
##'   ensemble members for each forcing (i.e., \eqn{n_i = } \code{nruns.X[i]}).
##' @param ctlruns an \eqn{m \times N} matrix of pre-industrial control run
##'   segments, assumed to have approximately zero mean.  Used to estimate the
##'   covariance structure of internal climate variability.
##' @param conf.level confidence level for marginal scaling factor intervals.
##'   Default is \code{0.90}.
##' @param sigma.X2 scalar inflation factor for the fingerprint pattern noise.
##'   If \code{NULL} (default), estimated automatically via \code{infFactor()}.
##' @param sigma.Z2 scalar inflation factor for the control run noise.
##'   If \code{NULL} (default), estimated automatically via \code{infFactor()}.
##' @param poly.grid an optional \eqn{M \times 3} matrix of candidate polynomial
##'   coefficients \eqn{(\ell_3, \ell_2, \ell_1)}, one per row.  Only
##'   monotonically admissible rows (i.e., \eqn{f} is positive and decreasing
##'   on \eqn{[0, \tau_{\max}]}) are evaluated.  If \code{NULL} (default), a
##'   \eqn{10 \times 10 \times 10} grid anchored to the Ledoit-Wolf reference
##'   value is generated and filtered automatically.
##'
##' @return A list with the following components:
##'   \itemize{
##'     \item \code{beta.optim}: length-\eqn{p} vector of scaling factor
##'       estimates at the optimal polynomial.
##'     \item \code{var.optim}: \eqn{p \times p} variance-covariance matrix
##'       at the optimal polynomial (\eqn{\hat{\Xi}(f_{\rm opt}) / N}).
##'     \item \code{ci.optim}: \eqn{2 \times p} matrix of marginal confidence
##'       interval bounds (rows: lower, upper).
##'     \item \code{sd.optim}: length-\eqn{p} vector of standard deviations at
##'       the optimal polynomial.
##'     \item \code{Xi.optim}: \eqn{p \times p} matrix \eqn{N \hat{\Xi}(f_{\rm opt})},
##'       the \eqn{N}-scaled asymptotic covariance.
##'     \item \code{sigma.X2}: scalar, estimated (or supplied) inflation factor
##'       for fingerprint pattern noise.
##'     \item \code{sigma.Z2}: scalar, estimated (or supplied) inflation factor
##'       for control run noise.
##'     \item \code{poly.optim}: length-3 vector \eqn{(\ell_3, \ell_2, \ell_1)}
##'       of the optimal polynomial coefficients.
##'     \item \code{trace_Xi}: length-\eqn{M} numeric vector of
##'       \eqn{\mathrm{tr}(\hat{\Xi}(f))} for every evaluated candidate.
##'       Useful for diagnosing the selection landscape.
##'     \item \code{z_scores}: length-\eqn{M} numeric vector of residual
##'       consistency z-scores for every evaluated candidate.
##'     \item \code{z_score.optim}: scalar, the residual consistency z-score at
##'       the optimal polynomial.  Values of \eqn{|z| > 2} suggest the
##'       polynomial family may not be flexible enough.
##'     \item \code{joint.cr}: a list giving the joint confidence region for
##'       \eqn{\beta}.  Contains \code{stat}
##'       (\eqn{\hat{\beta}^\top \hat{\Xi}^{-1}(f_{\rm opt})\hat{\beta}}),
##'       \code{bound} (\eqn{\chi^2_p(1-\alpha)} critical value), and
##'       \code{detected} (\code{TRUE} if the zero vector lies outside the
##'       joint confidence region).
##'     \item \code{poly.grid}: the filtered \eqn{M \times 3} matrix of
##'       admissible polynomial candidates that were evaluated.
##'   }
##'
##' @author Haoran Li, Yan Li
##' @keywords climate detection attribution fingerprinting regularization
##'   nonlinear shrinkage
##'
##' @references
##' \itemize{
##'   \item Li, H., Li, Y. (2025). Adaptable Fingerprinting with Nonlinear
##'     Shrinkage for Climate Change Detection and Attribution under Variance
##'     Heterogeneity.
##'   \item Li, Y., Chen, K., Yan, J. (2025). Regularized Fingerprinting with
##'     Linearly Optimal Weight Matrix in Detection and Attribution of Climate
##'     Change. \emph{Journal of Climate}.
##'   \item Ledoit, O., Wolf, M. (2004). A well-conditioned estimator for
##'     large-dimensional covariance matrices. \emph{Journal of Multivariate
##'     Analysis}, 88(2), 365--411.
##' }
##'
##' @examples
##' \dontrun{
##' data(globalDat)
##' ## Wrap the averaged fingerprint in a list of 1-row matrices to use
##' ## the simplified sigma.X2 = sigma.Z2 = 1 path.
##' X.list <- lapply(seq_len(ncol(globalDat$Xtilde)), function(i)
##'   matrix(globalDat$Xtilde[, i], nrow = 1))
##' result <- fingerprintNLS(
##'   X.tilde.list = X.list,
##'   Y            = globalDat$Y,
##'   nruns.X      = globalDat$nruns.X,
##'   ctlruns      = globalDat$ctlruns,
##'   sigma.X2     = 1,
##'   sigma.Z2     = 1
##' )
##' result$beta.optim
##' result$ci.optim
##' result$z_score.optim   ## residual consistency check
##' }
##'
##' @export fingerprintNLS
##' @importFrom stats cov qchisq qnorm
fingerprintNLS <- function(X.tilde.list, Y, nruns.X, ctlruns,
                           conf.level = 0.90,
                           sigma.X2 = NULL,
                           sigma.Z2 = NULL,
                           poly.grid = NULL) {
  Y <- as.matrix(Y)
  p <- length(X.tilde.list)
  N <- nrow(Y)
  ## compute ensemble means: N x p fingerprint matrix
  X.tilde <- sapply(X.tilde.list, colMeans)
  X.tilde <- as.matrix(X.tilde)
  if (is.null(colnames(X.tilde))) {
    ## inherit from list names if available (e.g. X.tilde.list = list(ANT=..., NAT=...))
    if (!is.null(names(X.tilde.list))) {
      colnames(X.tilde) <- names(X.tilde.list)
    } else {
      colnames(X.tilde) <- paste0("forcing ", seq_len(p))
    }
  }
  ## estimate inflation factors if not supplied
  if (is.null(sigma.X2) || is.null(sigma.Z2)) {
    sigma.hat <- infFactor(Y = Y, X.tilde.list = X.tilde.list,
                           Z = ctlruns, n.X = nruns.X)
    if (is.null(sigma.X2)) sigma.X2 <- sigma.hat["Sx"]
    if (is.null(sigma.Z2)) sigma.Z2 <- sigma.hat["Sz"]
  }
  ## build default polynomial candidate grid if not supplied
  if (is.null(poly.grid)) {
    ## anchor the l1 range to the Ledoit-Wolf reference value
    ## (requires lwRegcov.lam() from Covest_addition.R / R/Covest.R)
    l0 <- lwRegcov.lam(ctlruns)
    l3.seq <- seq(-2, 5, length = 10)
    l2.seq <- seq(-2, 5, length = 10)
    l1.seq <- seq(l0 * 0.5, l0 * 2, length = 10)
    poly.grid <- as.matrix(expand.grid(l3 = l3.seq, l2 = l2.seq, l1 = l1.seq))
    ## always include the pure linear shrinkage reference polynomial
    poly.grid <- rbind(poly.grid, c(1e-10, 1e-10, l0))
    ## keep only polynomials that are monotonically decreasing on [0, tau_max]
    lam.max <- max(eigen(cov(ctlruns), only.values = TRUE)$values) + 1e-3
    keep <- monotonic_increasing(poly.grid, lam.max)
    poly.grid <- poly.grid[keep, , drop = FALSE]
    if (nrow(poly.grid) == 0) {
      stop("No admissible polynomial candidates found; supply poly.grid manually.")
    }
  }
  ## guard: get_roots_b() divides by l3, so l3 = 0 produces NaN throughout
  if (any(poly.grid[, 1] == 0)) {
    stop("poly.grid contains rows with l3 = 0. The NLS method requires a ",
         "proper cubic polynomial (l3 != 0 in every candidate). Replace ",
         "zero entries with a small positive value such as 1e-8.")
  }
  ## evaluate all polynomial candidates: returns beta, Xi, trace_Xi, z_scores
  res <- Xi_polys(Y = Y, X.tilde = X.tilde, Z = ctlruns,
                  sigma.X2 = sigma.X2, sigma.Z2 = sigma.Z2,
                  n.X = nruns.X, poly.mat = poly.grid)
  ## select the optimal polynomial: minimises tr(Xi) = total asymptotic variance
  ind.opt <- which.min(res$trace.Xi)
  beta.opt <- res$beta[ind.opt, ]
  Xi.opt   <- res$Xi[ind.opt, , ]
  var.opt  <- Xi.opt / N
  ## marginal confidence intervals from asymptotic normality
  Z.crt  <- qnorm((1 - conf.level) / 2, lower.tail = FALSE)
  sd.opt <- sqrt(diag(var.opt))
  ci.opt <- rbind(beta.opt - Z.crt * sd.opt,
                  beta.opt + Z.crt * sd.opt)
  ## attach names
  names(beta.opt) <- names(sd.opt) <- colnames(X.tilde)
  rownames(ci.opt) <- c("lower", "upper")
  colnames(ci.opt) <- colnames(X.tilde)
  ## joint confidence region
  jcr.stat  <- as.numeric(t(beta.opt) %*% solve(var.opt) %*% beta.opt)
  jcr.bound <- qchisq(conf.level, df = p)
  list(
    beta.optim    = beta.opt,
    var.optim     = var.opt,
    ci.optim      = ci.opt,
    sd.optim      = sd.opt,
    Xi.optim      = Xi.opt,
    sigma.X2      = as.numeric(sigma.X2),
    sigma.Z2      = as.numeric(sigma.Z2),
    poly.optim    = poly.grid[ind.opt, ],
    trace_Xi      = res$trace.Xi,
    z_scores      = res$z.scores,
    z_score.optim = res$z.scores[ind.opt],
    joint.cr      = list(
      stat     = jcr.stat,
      bound    = jcr.bound,
      detected = jcr.stat >= jcr.bound
    ),
    poly.grid     = poly.grid
  )
}


## Estimate variance inflation factors sigma.X2 and sigma.Z2.
##
## Performs a grid search over candidate values of sigma.X2 (the measurement
## error inflation factor) by minimising the discrepancy between the observed
## within-ensemble variance and a model-implied quantity.  sigma.Z2 is then
## derived as sigma.X2 * Vz / Vx where Vz is the mean diagonal of the sample
## covariance of the control runs.
##
## @param Y       N x 1 response matrix.
## @param X.tilde.list  list of p matrices, each nᵢ x N (raw ensemble members).
## @param Z       m x N control run matrix.
## @param n.X     length-p vector of ensemble sizes.
## @param W       optional N x N weight matrix.  If NULL, uses lwRegcov(cov(Z)).
## @param seq.range  length-2 vector: search range for a (default c(0.1, 5)).
## @param seq.length  integer: number of candidate a values to try (default 100).
## @return named numeric vector c(Sx = sigma.X2, Sz = sigma.Z2).
infFactor <- function(Y, X.tilde.list, Z, n.X,
                      W = NULL,
                      seq.range = c(0.1, 5),
                      seq.length = 100) {
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  N   <- nrow(Y)
  m   <- nrow(Z)
  p   <- length(X.tilde.list)
  D.mat  <- diag(1 / n.X)
  D.sinv <- diag(sqrt(n.X))
  S <- cov(Z)
  ## ensemble means: N x p matrix
  X.tilde <- sapply(X.tilde.list, colMeans)
  ## Vx: pooled within-ensemble variance across all forcings and observations
  VX <- sum(sapply(X.tilde.list, function(mat) {
    sum(mat^2) - nrow(mat) * sum(colMeans(mat)^2)
  })) / N / (sum(n.X) - p)
  ## Vz: mean diagonal of the sample covariance of the control runs
  VZ <- mean(diag(S))
  ## default weight matrix: Ledoit-Wolf linear shrinkage
  if (is.null(W)) {
    W <- lwRegcov(S)
  }
  W.inv <- solve(W)
  ## precomputed cross-products
  XWX <- t(X.tilde) %*% W.inv %*% X.tilde
  XWY <- t(X.tilde) %*% W.inv %*% Y
  YWY <- t(Y)       %*% W.inv %*% Y
  ## grid search over candidate values of a = sigma.X2
  cand.a <- seq(seq.range[1], seq.range[2], length = seq.length)
  tmp <- vapply(cand.a, function(a) {
    ## augmented matrix at candidate a
    MtM <- rbind(
      cbind(D.sinv %*% XWX %*% D.sinv,
            -sqrt(a) * D.sinv %*% XWY),
      cbind(-sqrt(a) * t(XWY) %*% D.sinv,
            a * YWY)
    )
    ## smallest eigenvector gives the TLS solution at this a
    ev    <- eigen(MtM)$vectors[, p + 1]
    beta  <- D.sinv %*% (ev[seq_len(p)] / ev[p + 1] / sqrt(a))
    ## objective: discrepancy between model-implied and observed within-ensemble var
    discr <- abs(a * sum((Y - X.tilde %*% beta)^2) / N /
                   (1 + a * as.numeric(t(beta) %*% D.mat %*% beta)) - VX)
    c(a = a, discr = discr)
  }, numeric(2))
  sigma.xs         <- tmp["a",    which.min(tmp["discr", ])]
  sigma.zs         <- sigma.xs * VZ / VX
  names(sigma.xs)  <- "Sx"
  names(sigma.zs)  <- "Sz"
  c(sigma.xs, sigma.zs)
}


## Compute beta, Xi, trace_Xi, and residual z-scores for all polynomial candidates.
##
## This is the main computational engine called by fingerprintNLS().  For each
## row of poly.mat, it evaluates the TLS estimator and the RMT-corrected
## asymptotic covariance matrix.
##
## @param Y        N x 1 response matrix.
## @param X.tilde  N x p matrix of ensemble-mean fingerprint patterns.
## @param Z        m x N control run matrix (assumed mean-zero).
## @param sigma.X2 scalar inflation factor for measurement error.
## @param sigma.Z2 scalar inflation factor for control run noise.
## @param n.X      length-p vector of ensemble sizes.
## @param poly.mat M x 3 matrix of polynomial coefficients (l3, l2, l1).
## @return list with beta (M x p), Xi (M x p x p array), trace.Xi (length-M),
##   z.scores (length-M), Delta (M x p x p), Omega (M x p x p), K (length-M).
Xi_polys <- function(Y, X.tilde, Z, sigma.X2, sigma.Z2, n.X, poly.mat) {
  N     <- nrow(Y)
  p     <- ncol(X.tilde)
  m     <- nrow(Z)
  npoly <- nrow(poly.mat)
  D.mat <- diag(1 / n.X)
  ## raw second-moment matrix of control runs (assumes zero mean)
  S <- (1 / m) * t(Z) %*% Z
  ## project X and Y onto eigenspace of S; inflate X by 1/sqrt(sigma.X2)
  tmp        <- eigen_S_proj_X(X.tilde / sqrt(sigma.X2), Y, S)
  eigenvalues <- tmp$eigenvalues
  Q.X         <- tmp$Q.X
  Q.Y         <- tmp$Q.Y
  QQ.list     <- tmp$QQ.list
  ## TLS estimates for all polynomial candidates
  tmp.beta     <- beta_poly(poly.mat, Q.X, Q.Y, n.X, N, m, eigenvalues)
  beta.list    <- tmp.beta$beta.list
  ## divide residuals by N to get mean residual per observation
  residual.list <- tmp.beta$residual.list / N
  ## RMT auxiliary quantities for all candidates
  tmp.dok  <- Delta_Omega_K_poly(poly.mat, QQ.list, n.X, N, m,
                                 eigenvalues, sigma.Z2 = sigma.Z2)
  Delta.list <- tmp.dok$Delta
  Omega.list <- tmp.dok$Omega
  K.list     <- tmp.dok$K
  PI.list    <- tmp.dok$PI
  ## preallocate outputs
  Xi.array <- array(NA_real_, dim = c(npoly, p, p))
  trace.Xi <- numeric(npoly)
  z.scores <- numeric(npoly)
  for (ss in seq_len(npoly)) {
    beta.f  <- beta.list[ss, ]
    Delta.f <- Delta.list[ss, , ]
    Omega.f <- Omega.list[ss, , ]
    K.f     <- K.list[ss]
    PI.f    <- PI.list[ss]
    inv.Delta <- solve(Delta.f)
    mid.mat   <- solve(diag(n.X) + beta.f %*% t(beta.f))
    scale1    <- as.numeric(1 + t(beta.f) %*% D.mat %*% beta.f)
    Xi.array[ss, , ] <- scale1 * inv.Delta %*%
      (Omega.f + K.f * mid.mat) %*% inv.Delta
    trace.Xi[ss] <- sum(diag(Xi.array[ss, , ]))
    ## residual consistency z-score
    exp.res <- scale1 * PI.f
    var.res <- 2 / N * scale1 * K.f
    z.scores[ss] <- (residual.list[ss] - exp.res) / sqrt(var.res)
  }
  ## back-adjust for inflation factor sigma.X2
  list(beta     = beta.list / sqrt(sigma.X2),
       Xi       = Xi.array  / sigma.X2,
       trace.Xi = trace.Xi  / sigma.X2,
       z.scores = z.scores,
       Delta    = Delta.list,
       Omega    = Omega.list,
       K        = K.list)
}


## Residual consistency test for a single polynomial shrinkage function.
##
## Computes a z-score testing whether the observed weighted residuals are
## consistent with their theoretical expectation under the assumed model.
## Under the null (model correct), z -> N(0,1).  |z| > 2 indicates possible
## model misspecification.
##
## @param Y        N x 1 response matrix.
## @param X.tilde  N x p fingerprint pattern matrix.
## @param Z        m x N control run matrix (assumed mean-zero).
## @param poly.coef  matrix with one or more rows, each giving (l3, l2, l1).
## @param sigma.X2  scalar inflation factor for measurement error.
## @param sigma.Z2  scalar inflation factor for control run noise.
## @param n.X      length-p vector of ensemble sizes.
## @return list with p.value (z-score vector) and mse.residual.list.
residual_test <- function(Y, X.tilde, Z, poly.coef,
                          sigma.X2, sigma.Z2, n.X) {
  N   <- nrow(Y)
  p   <- ncol(X.tilde)
  m   <- nrow(Z)
  D.mat <- diag(1 / n.X)
  if (!is.matrix(poly.coef)) poly.coef <- matrix(poly.coef, nrow = 1)
  S <- (1 / m) * t(Z) %*% Z
  tmp         <- eigen_S_proj_X(X.tilde / sqrt(sigma.X2), Y, S)
  eigenvalues <- tmp$eigenvalues
  Q.X         <- tmp$Q.X
  Q.Y         <- tmp$Q.Y
  QQ.list     <- tmp$QQ.list
  z.score.list     <- numeric(nrow(poly.coef))
  mse.residual.list <- numeric(nrow(poly.coef))
  for (i in seq_len(nrow(poly.coef))) {
    coef.i <- poly.coef[i, , drop = FALSE]
    beta.i <- beta_poly(coef.i, Q.X, Q.Y, n.X, N, m, eigenvalues)$beta.list
    ## polynomial weights applied to eigenvalues
    poly.wt <- 1 / (coef.i[1] * eigenvalues^3 +
                      coef.i[2] * eigenvalues^2 +
                      coef.i[3] * eigenvalues  + 1)
    tmp.cov <- Delta_Omega_K(coef.i, QQ.list, n.X, N, m,
                             eigenvalues, sigma.Z2 = sigma.Z2)
    Delta.i <- tmp.cov$Delta
    PI.i    <- tmp.cov$PI
    Omega.i <- tmp.cov$Omega
    K.i     <- tmp.cov$K
    mse.res  <- mean(poly.wt * (Q.Y - Q.X %*% as.vector(beta.i))^2)
    scale1   <- as.numeric(1 + beta.i %*% D.mat %*% t(beta.i))
    exp.res  <- scale1 * PI.i
    var.res  <- 2 / N * scale1 * K.i
    z.score.list[i]      <- (mse.res - exp.res) / sqrt(var.res)
    mse.residual.list[i] <- mse.res
  }
  list(p.value          = z.score.list,
       mse.residual.list = mse.residual.list)
}


## TLS estimate given an explicit weight matrix.
##
## Utility function for users who wish to supply their own weight matrix
## rather than using the polynomial-family search.
##
## @param Y        N x 1 response matrix.
## @param X.tilde  N x p fingerprint pattern matrix.
## @param Z        m x N control run matrix (unused here, kept for interface
##   consistency with other internal functions).
## @param sigma.X2 scalar inflation factor for measurement error.
## @param n.X      length-p vector of ensemble sizes.
## @param weight   N x N positive definite weight matrix.
## @return length-p vector of scaling factor estimates.
beta_est <- function(Y, X.tilde, Z, sigma.X2, n.X, weight) {
  N      <- length(Y)
  p      <- ncol(X.tilde)
  W.inv  <- solve(weight)
  ## scale X columns by sqrt(n.X) and apply inflation adjustment
  scaled.X <- X.tilde * t(sqrt(n.X) %*% matrix(1, 1, N))
  bind.XY  <- cbind(scaled.X, sqrt(sigma.X2) * Y)
  MtM      <- t(bind.XY) %*% W.inv %*% bind.XY
  lam.min  <- tail(eigen(MtM)$values, 1)
  if (length(n.X) != 1) {
    Dn.X <- diag(sqrt(n.X))
  } else {
    Dn.X <- sqrt(n.X)
  }
  as.vector(Dn.X %*%
    as.vector(solve(MtM[seq_len(p), seq_len(p)] - lam.min * diag(p)) %*%
                MtM[p + 1, seq_len(p)])) / sqrt(sigma.X2)
}


## TLS estimates for a family of polynomial shrinkage functions.
##
## For each row of poly.mat, computes the TLS estimator and the residual
## statistic (1 + ||beta||^2) * lambda_min.  Works entirely in the
## eigenspace of S for efficiency.
##
## @param poly.mat   M x 3 matrix of coefficients (l3, l2, l1).
## @param Q.X        N x p rotated fingerprint matrix t(U) %*% X.
## @param Q.Y        N x 1 rotated response vector t(U) %*% Y.
## @param n.X        length-p ensemble sizes.
## @param N          integer, observation dimension.
## @param m          integer, number of control run segments.
## @param eigenvalues  length-N eigenvalues of S.
## @return list with beta.list (M x p matrix) and residual.list (length-M).
beta_poly <- function(poly.mat, Q.X, Q.Y, n.X, N, m, eigenvalues) {
  npoly     <- nrow(poly.mat)
  p         <- ncol(Q.X)
  scaled.X  <- Q.X * t(sqrt(n.X) %*% matrix(1, 1, N))
  if (length(n.X) != 1) {
    Dn.X <- diag(sqrt(n.X))
  } else {
    Dn.X <- sqrt(n.X)
  }
  bind.XY <- cbind(scaled.X, Q.Y)
  ## outer products of each row of [scaled.X | Q.Y], used to build MtM
  YY.list <- lapply(seq_len(nrow(bind.XY)),
                    function(i) bind.XY[i, ] %*% t(bind.XY[i, ]))
  ## for each polynomial, compute weighted MtM and TLS solution
  result <- apply(poly.mat, 1, function(x) {
    ## polynomial weights: f(tau_i) for each eigenvalue
    poly.wt <- 1 / (x[1] * eigenvalues^3 +
                      x[2] * eigenvalues^2 +
                      x[3] * eigenvalues  + 1)
    MtM     <- Reduce("+", Map("*", poly.wt, YY.list))
    lam.min <- tail(eigen(MtM)$values, 1)
    beta    <- as.vector(
      solve(MtM[seq_len(p), seq_len(p)] - lam.min * diag(p)) %*%
        MtM[p + 1, seq_len(p)]
    )
    res <- (sum(beta^2) + 1) * lam.min
    c(beta, res)
  })
  residual.list <- result[p + 1, ]
  ## back-transform to original scale
  beta.list <- t(Dn.X %*% result[seq_len(p), ])
  list(beta.list = beta.list, residual.list = residual.list)
}


## Compute Delta, Omega, K, and PI for all polynomial candidates (batch).
##
## Thin wrapper that loops over rows of poly.mat and calls Delta_Omega_K().
##
## @param poly.mat    M x 3 coefficient matrix.
## @param QQ.list     list of N outer-product matrices (p x p each).
## @param n.X         length-p ensemble sizes.
## @param N           observation dimension.
## @param m           number of control run segments.
## @param eigenvalues  length-N eigenvalues of S.
## @param sigma.Z2    scalar inflation factor for control run noise.
## @return list with Delta (M x p x p), Omega (M x p x p), K (length-M),
##   PI (length-M).
Delta_Omega_K_poly <- function(poly.mat, QQ.list, n.X, N, m,
                                eigenvalues, sigma.Z2) {
  npoly  <- nrow(poly.mat)
  p      <- ncol(QQ.list[[1]])
  Delta.arr <- array(NA_real_, dim = c(npoly, p, p))
  Omega.arr <- array(NA_real_, dim = c(npoly, p, p))
  K.vec     <- numeric(npoly)
  PI.vec    <- numeric(npoly)
  for (ii in seq_len(npoly)) {
    tmp <- Delta_Omega_K(poly.mat[ii, ], QQ.list, n.X, N, m,
                         eigenvalues, sigma.Z2 = sigma.Z2)
    Delta.arr[ii, , ] <- tmp$Delta
    Omega.arr[ii, , ] <- tmp$Omega
    K.vec[ii]          <- tmp$K
    PI.vec[ii]         <- tmp$PI
  }
  list(Delta = Delta.arr, Omega = Omega.arr, K = K.vec, PI = PI.vec)
}


## Compute Delta(f), Omega(f), K(f), and PI(f) for a single polynomial.
##
## Uses the partial fraction decomposition of f to express the RMT auxiliary
## quantities as sums over the three roots of the denominator polynomial.
## All complex quantities are reduced to their real parts before returning.
##
## @param coef       length-3 vector (l3, l2, l1) or 1 x 3 matrix.
## @param QQ.list    list of N outer-product matrices (p x p each).
## @param n.X        length-p ensemble sizes (passed as D = 1/n.X to delta_hat).
## @param N          observation dimension.
## @param m          number of control run segments.
## @param eigenvalues  length-N eigenvalues of S.
## @param sigma.Z2   scalar inflation factor for control run noise.
## @return list with Delta (p x p), Omega (p x p), K (scalar), PI (scalar).
Delta_Omega_K <- function(coef, QQ.list, n.X, N, m, eigenvalues, sigma.Z2) {
  coef    <- as.vector(coef)
  tmp.pfr <- get_roots_b(coef)
  r <- tmp.pfr$r
  b <- tmp.pfr$b
  k <- length(r)
  D <- 1 / n.X
  ## theta and theta' for each root
  tmp.theta     <- theta_hat(r, N, m, eigenvalues)
  theta.vec     <- tmp.theta[1, ]
  theta.prime.vec <- tmp.theta[2, ]
  ## delta and delta' matrices for each root
  tmp.delta     <- delta_hat(r, QQ.list, D, N, m, eigenvalues,
                             theta.vec, theta.prime.vec, sigma.Z2)
  delta.list        <- tmp.delta$delta
  delta.prime.list  <- tmp.delta$delta.prime
  ## Delta(f) = sum_j b_j * delta_j
  Delta.f <- Reduce("+", Map("*", b, delta.list))
  ## PI(f) = sum_j b_j * (m/N) * (theta_j - 1)
  PI.f    <- Reduce("+", Map("*", b, (m / N) * (theta.vec - 1)))
  ## K(f) = sum_i sum_{j>=i} b_i b_j * kappa(i,j) * (1 + (i!=j))
  K.f <- 0
  for (i in seq_len(k)) {
    for (j in i:k) {
      kap <- kappa_hat(i, j, r, theta.vec, theta.prime.vec, N, m)
      K.f <- K.f + kap * b[i] * b[j] * ifelse(i == j, 1, 2)
    }
  }
  ## Omega(f) = sum_i sum_{j>=i} b_i b_j * omega(i,j) * (1 + (i!=j))
  Omega.f <- 0
  for (i in seq_len(k)) {
    for (j in i:k) {
      omg <- omega_hat(i, j, r, theta.vec, delta.list, delta.prime.list)
      Omega.f <- Omega.f + omg * b[i] * b[j] * ifelse(i == j, 1, 2)
    }
  }
  list(Delta = Re(Delta.f),
       PI    = Re(PI.f)    / sigma.Z2,
       Omega = Re(Omega.f) / sigma.Z2,
       K     = Re(K.f)     / sigma.Z2^2)
}


## Stieltjes transform of the empirical spectral distribution: mean(1/(tau - z)).
##
## Evaluated at a complex number z (below the real axis) representing a root
## of the polynomial denominator.
##
## @param z         complex scalar (a root of the denominator polynomial).
## @param emp.eig   length-N vector of sample eigenvalues of S.
## @return complex scalar.
m_hat <- function(z, emp.eig) {
  mean(1 / (emp.eig - z))
}


## Derivative of the Stieltjes transform: mean(1/(tau - z)^2).
##
## @param z         complex scalar.
## @param emp.eig   length-N vector of sample eigenvalues of S.
## @return complex scalar.
m_hat_prime <- function(z, emp.eig) {
  mean(1 / (emp.eig - z)^2)
}


## Compute theta_hat and its derivative for each root of the polynomial.
##
## theta_hat(z) = 1 / (1 - gamma - gamma * z * m_hat(z)) where gamma = N/m.
## These are RMT correction scalars derived from the Marchenko-Pastur law.
##
## @param r           complex vector of polynomial roots.
## @param N           observation dimension.
## @param m           number of control run segments.
## @param eigenvalues  length-N eigenvalues of S.
## @return 2 x length(r) matrix: row 1 = theta, row 2 = theta'.
theta_hat <- function(r, N, m, eigenvalues) {
  gamma <- N / m
  sapply(r, function(z) {
    m.z      <- m_hat(z, eigenvalues)
    m.prime  <- m_hat_prime(z, eigenvalues)
    theta    <- 1 / (1 - gamma - gamma * z * m.z)
    theta.pr <- gamma * (m.z + z * m.prime) * theta^2
    c(theta = theta, theta.prime = theta.pr)
  })
}


## Compute delta_hat and its derivative matrices for each polynomial root.
##
## delta_hat(z) = (1/N) sum_i [1/(tau_i - z)] Q_X[i,] Q_X[i,]^T
##                - (m/N)(theta(z) - 1) D / sigma.Z2
##
## @param r                 complex vector of polynomial roots.
## @param QQ.list           list of N outer-product p x p matrices.
## @param D                 length-p vector 1/n.X (converted to diag matrix inside).
## @param N                 observation dimension.
## @param m                 number of control run segments.
## @param eigenvalues        length-N eigenvalues of S.
## @param theta.vec         length-k vector of theta values at each root.
## @param theta.prime.vec   length-k vector of theta' values at each root.
## @param sigma.Z2          scalar inflation factor for control run noise.
## @return list with delta (list of k complex p x p matrices) and
##   delta.prime (list of k complex p x p matrices).
delta_hat <- function(r, QQ.list, D, N, m, eigenvalues,
                      theta.vec, theta.prime.vec, sigma.Z2) {
  if (is.null(dim(D))) D <- diag(D)
  delta.list <- lapply(seq_along(r), function(i) {
    z     <- r[i]
    ss    <- 1 / (eigenvalues - z)
    term1 <- (1 / N) * Reduce("+", Map("*", ss, QQ.list))
    term2 <- -(m / N) * (theta.vec[i] - 1) * D / sigma.Z2
    term1 + term2
  })
  delta.prime.list <- lapply(seq_along(r), function(i) {
    z      <- r[i]
    ss.sq  <- 1 / (eigenvalues - z)^2
    term3  <- (1 / N) * Reduce("+", Map("*", ss.sq, QQ.list))
    term4  <- -(m / N) * theta.prime.vec[i] * D / sigma.Z2
    term3 + term4
  })
  list(delta = delta.list, delta.prime = delta.prime.list)
}


## Compute kappa_hat(i, j): scalar cross-term for K(f).
##
## kappa_hat(i,i) = (1/gamma) * theta_i^2 * (theta_i + r_i * theta'_i - 1)
## kappa_hat(i,j) = (1/gamma) * theta_i * theta_j *
##                  ((r_i theta_i - r_j theta_j) / (r_i - r_j) - 1)  for i != j
##
## @param i,j              root indices (integers).
## @param r                complex vector of polynomial roots.
## @param theta.vec        length-k theta values.
## @param theta.prime.vec  length-k theta' values.
## @param N                observation dimension.
## @param m                number of control run segments.
## @return complex scalar.
kappa_hat <- function(i, j, r, theta.vec, theta.prime.vec, N, m) {
  gamma <- N / m
  if (i == j) {
    (1 / gamma) * theta.vec[i]^2 *
      (theta.vec[i] + r[i] * theta.prime.vec[i] - 1)
  } else {
    (1 / gamma) * theta.vec[i] * theta.vec[j] *
      ((r[i] * theta.vec[i] - r[j] * theta.vec[j]) /
         (r[i] - r[j]) - 1)
  }
}


## Compute omega_hat(i, j): p x p matrix cross-term for Omega(f).
##
## omega_hat(i,i) = theta_i^2 * (delta_i + r_i * delta'_i)
## omega_hat(i,j) = theta_i * theta_j *
##                  (r_i * delta_i - r_j * delta_j) / (r_i - r_j)  for i != j
##
## @param i,j               root indices (integers).
## @param r                 complex vector of polynomial roots.
## @param theta.vec         length-k theta values.
## @param delta.list        list of k complex p x p delta matrices.
## @param delta.prime.list  list of k complex p x p delta' matrices.
## @return complex p x p matrix.
omega_hat <- function(i, j, r, theta.vec, delta.list, delta.prime.list) {
  if (i == j) {
    theta.vec[i]^2 * (delta.list[[i]] + r[i] * delta.prime.list[[i]])
  } else {
    theta.vec[i] * theta.vec[j] *
      (r[i] * delta.list[[i]] - r[j] * delta.list[[j]]) / (r[i] - r[j])
  }
}


## Eigendecompose S and project X.tilde and Y onto the eigenbasis.
##
## Computes S = U Lambda U^T, then Q_X = U^T X.tilde and Q_Y = U^T Y.
## Also builds QQ.list: a list of N outer-product matrices Q_X[i,] Q_X[i,]^T.
##
## @param X.tilde   N x p fingerprint pattern matrix (may be pre-scaled by
##   1/sqrt(sigma.X2)).
## @param Y         N x 1 response vector.
## @param S         N x N covariance matrix to eigendecompose.
## @return list with eigenvalues (length-N), Q.X (N x p), Q.Y (N x 1),
##   QQ.list (list of N p x p matrices).
eigen_S_proj_X <- function(X.tilde, Y, S) {
  tmp         <- eigen(S)
  eigenvalues <- tmp$values
  Q.X         <- t(tmp$vectors) %*% X.tilde
  Q.Y         <- t(tmp$vectors) %*% Y
  QQ.list     <- lapply(seq_len(nrow(X.tilde)),
                        function(i) Q.X[i, ] %*% t(Q.X[i, ]))
  list(eigenvalues = eigenvalues, Q.X = Q.X, Q.Y = Q.Y, QQ.list = QQ.list)
}


## Partial fraction decomposition of a cubic polynomial denominator.
##
## Given f(x) = 1 / (l3*x^3 + l2*x^2 + l1*x + 1), finds the roots
## r1, r2, r3 and residues b1, b2, b3 such that
## f(x) = b1/(x - r1) + b2/(x - r2) + b3/(x - r3).
##
## @param coef  length-3 vector (l3, l2, l1).
## @return list with r (length-3 complex roots) and b (length-3 complex
##   residues).
get_roots_b <- function(coef) {
  ## polyroot expects coefficients from constant to highest degree
  r  <- polyroot(c(1, coef[length(coef):1]))
  l3 <- coef[1]
  r1 <- r[1]; r2 <- r[2]; r3 <- r[3]
  b1 <- 1 / (l3 * (r1 - r2) * (r1 - r3))
  b2 <- 1 / (l3 * (r2 - r1) * (r2 - r3))
  b3 <- 1 / (l3 * (r3 - r1) * (r3 - r2))
  list(b = c(b1, b2, b3), r = c(r1, r2, r3))
}


## Screen a matrix of cubic polynomial candidates for monotone admissibility.
##
## A candidate f(x) = 1/(l3*x^3 + l2*x^2 + l1*x + 1) is admissible if the
## denominator polynomial is positive and monotonically increasing on [0, alpha],
## ensuring that f itself is positive and decreasing (the "more variability ->
## less weight" property required for a valid precision matrix).
##
## The minimum of the derivative l'(x) = 3*l3*x^2 + 2*l2*x + l1 on [0, alpha]
## is checked at the three candidates x = 0, x = alpha, and the interior vertex
## x* = -l2/(3*l3) if it lies in the interval.
##
## @param A      M x 3 matrix of coefficients (l3, l2, l1), one candidate per row.
## @param alpha  positive scalar: upper end of the interval (typically tau_max).
## @return logical vector of length M: TRUE if the polynomial is admissible.
monotonic_increasing <- function(A, alpha) {
  ## derivative at x=0: l'(0) = l1
  cond1 <- A[, 3] > 0
  ## derivative at x=alpha
  cond2 <- (3 * A[, 1] * alpha^2 + 2 * A[, 2] * alpha + A[, 3]) > 0
  ## interior vertex x* = -l2/(3*l3); if in (0, alpha), check derivative there
  A.safe      <- A
  A.safe[A.safe[, 1] == 0, 1] <- 1e-10
  x.vertex    <- -A.safe[, 2] / (3 * A.safe[, 1])
  x.check     <- ifelse(x.vertex > 0 & x.vertex < alpha, x.vertex, 0)
  cond3 <- (3 * A[, 1] * x.check^2 + 2 * A[, 2] * x.check + A[, 3]) > 0
  cond1 & cond2 & cond3
}

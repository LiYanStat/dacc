# Internal helpers for polynomial fingerprinting.
# Not exported. Numerical inputs are passed explicitly by fingerprintPoly().

# Validate a finite integer scalar.
.fp_integer_scalar <- function(x, lower) is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x >= lower && x == floor(x)

# Validate a finite positive scalar.
.fp_positive <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x > 0

# Validate a finite numeric matrix.
.fp_finite_matrix <- function(x) is.matrix(x) && is.numeric(x) && all(is.finite(x))

# Evaluate a polynomial in ascending coefficient order (Horner scheme).
.fp_poly <- function(a, x) {
    out <- rep(0, length(x))
    for (j in rev(seq_along(a))) out <- out * x + a[j]
    out
}

# Check monotonicity, boundedness and root separation.
.fp_admissible <- function(a, degree, psi1, psi2, psi3) {
    if (a[length(a)] == 0)
        return(FALSE)
    da <- seq_len(degree) * a[-1]
    critical <- c(0, psi1)
    if (degree > 1) {
        rr <- polyroot(seq_len(degree - 1) * da[-1])
        critical <- c(critical, Re(rr[abs(Im(rr)) < 1e-08 & Re(rr) > 0 & Re(rr) <
            psi1]))
    }
    if (min(.fp_poly(da, critical)) < -1e-10 || .fp_poly(a, psi1) > psi2)
        return(FALSE)
    roots <- polyroot(a)
    if (any(!is.finite(roots)))
        return(FALSE)
    if (degree > 1 && min(Mod(outer(roots, roots, "-"))[upper.tri(matrix(0, degree,
        degree))]) < psi3)
        return(FALSE)
    TRUE
}

# Compute weighted TLS in the sample covariance eigenbasis.
.fp_tls <- function(a, sx, ev, Xt, yt, N, mruns, D) {
    w <- .fp_poly(a, ev)
    G <- crossprod(Xt, Xt/w)/N
    b <- crossprod(Xt, yt/w)/N
    A <- rbind(cbind(G * outer(sqrt(mruns), sqrt(mruns)), sqrt(sx * mruns) * b),
        c(sqrt(sx * mruns) * b, sx * sum(yt^2/w)/N))
    lam <- min(eigen(A, symmetric = TRUE, only.values = TRUE)$values)
    beta <- as.vector(solve(G - lam * D, b))
    list(beta = beta, lambda = lam, w = w)
}

# Discard negligible imaginary roundoff after checking stability.
.fp_real <- function(x) {
    if (any(!is.finite(x)) || max(abs(Im(x))) > 1e-06 * max(1, max(abs(Re(x)))))
        stop("Unstable partial fractions")
    Re(x)
}

# Estimate asymptotic covariance and trace terms (Appendix A).
.fp_moments <- function(a, sx, sz, fit, degree, N, m, ev, Xt, D, p, mruns) {
    roots <- polyroot(a)
    weights <- 1/.fp_poly(seq_len(degree) * a[-1], roots)
    rho <- N/m
    objs <- lapply(roots, function(r) {
        h <- 1/(ev - r)
        mu <- mean(h)
        mup <- mean(h^2)
        theta <- 1/(1 - rho - rho * r * mu)
        thetap <- rho * (mu + r * mup) * theta^2
        pi <- (theta - 1)/rho
        pip <- thetap/rho
        delta <- crossprod(Xt, Xt * h)/N - sx/sz * pi * D
        deltap <- crossprod(Xt, Xt * h^2)/N - sx/sz * pip * D
        list(t = theta, pi = pi, pip = pip, d = delta, dp = deltap)
    })
    Delta <- matrix(0 + (0+0i), p, p)
    Omega <- Delta
    K <- 0 + (0+0i)
    Pi <- 0 + (0+0i)
    for (i in seq_len(degree)) {
        u <- objs[[i]]
        Delta <- Delta + weights[i] * u$d
        Pi <- Pi + weights[i] * u$pi/sz
        for (j in seq_len(degree)) {
            v <- objs[[j]]
            if (i == j) {
                kap <- u$t^2 * (roots[i] * u$pip + u$pi)
                om <- u$t^2 * (u$d + roots[i] * u$dp)
            }
            else {
                kap <- u$t * v$t * (roots[i] * u$pi - roots[j] * v$pi)/(roots[i] -
                  roots[j])
                om <- u$t * v$t * (roots[i] * u$d - roots[j] * v$d)/(roots[i] -
                  roots[j])
            }
            K <- K + weights[i] * weights[j] * kap/sz^2
            Omega <- Omega + weights[i] * weights[j] * om/sz
        }
    }
    Delta <- .fp_real(Delta)
    Omega <- .fp_real(Omega)
    K <- .fp_real(K)
    Pi <- .fp_real(Pi)
    chol(Delta)
    chol(Omega)
    if (K <= 0 || Pi <= 0)
        stop("Nonpositive trace estimates")
    inv <- solve(Delta)
    beta <- fit$beta
    scale <- 1 + sx * sum(beta^2/mruns)
    Xi <- scale * inv %*% (Omega + sx * K * solve(diag(mruns, p) + sx * tcrossprod(beta))) %*%
        inv
    Xi <- (Xi + t(Xi))/2
    chol(Xi)
    list(Xi = Xi, Pi = Pi, K = K, scale = scale)
}

# Evaluate the variance inflation matching loss.
.fp_loss <- function(s, a, ev, Xt, yt, N, mruns, D, Y, Xtilde, VX) {
    fit <- tryCatch(.fp_tls(a, s, ev = ev, Xt = Xt, yt = yt, N = N, mruns = mruns,
        D = D), error = function(e) NULL)
    if (is.null(fit))
        return(Inf)
    abs(s * sum((Y - Xtilde %*% fit$beta)^2)/(N * (1 + s * sum(fit$beta^2/mruns))) -
        VX)
}

# Search the variance inflation grid and recover both coefficients.
.fp_estimate <- function(a, ev, Xt, yt, N, mruns, D, inflation.grid, Y, Xtilde, VX,
    VZ) {
    grid <- inflation.grid
    if (is.null(grid)) {
        w <- .fp_poly(a, ev)
        U <- sweep(Xt, 2, sqrt(mruns), "*")
        A0 <- crossprod(U, U/w)/N
        b <- crossprod(U, yt/w)/N
        cc <- sum(yt^2/w)/N
        limit <- min(eigen(A0 - tcrossprod(b)/cc, symmetric = TRUE, only.values = TRUE)$values)
        if (!is.finite(limit) || limit <= 0)
            stop("Automatic factor-grid bound requires linearly independent columns of [Xtilde, Y]; this is a design condition, not a positive-definiteness requirement on the control covariance")
        upper <- 1
        while (.fp_tls(a, upper, ev = ev, Xt = Xt, yt = yt, N = N, mruns = mruns,
            D = D)$lambda < 0.95 * limit && upper < 2^20) upper <- upper * 2
        for (attempt in seq_len(9)) {
            grid <- sort(unique(c(seq(0, upper, length.out = 201), exp(seq(log(upper) -
                log(1e+08), log(upper), length.out = 201)))))
            losses <- vapply(grid, .fp_loss, numeric(1), a = a, ev = ev, Xt = Xt,
                yt = yt, N = N, mruns = mruns, D = D, Y = Y, Xtilde = Xtilde, VX = VX)
            idx <- which.min(losses)
            if (idx < 0.95 * length(grid))
                break
            upper <- upper * 2
        }
    }
    losses <- vapply(grid, .fp_loss, numeric(1), a = a, ev = ev, Xt = Xt, yt = yt,
        N = N, mruns = mruns, D = D, Y = Y, Xtilde = Xtilde, VX = VX)
    idx <- which.min(losses)
    if (!is.finite(losses[idx]) || grid[idx] <= 0)
        stop("No positive, finite inflation-coefficient estimate on grid")
    if (idx >= 0.95 * length(grid))
        warning("Inflation estimate near upper grid endpoint; enlarge inflation.grid")
    c(x = grid[idx], z = grid[idx] * VZ/VX)
}

# Select the admissible polynomial with smallest total variance.
.fp_select <- function(inflation, candidates, ev, Xt, yt, N, mruns, D, degree, m,
    p) {
    best <- NULL
    score <- Inf
    for (a in candidates) {
        obj <- tryCatch({
            fit <- .fp_tls(a, inflation[1], ev = ev, Xt = Xt, yt = yt, N = N, mruns = mruns,
                D = D)
            mom <- .fp_moments(a, inflation[1], inflation[2], fit, degree = degree,
                N = N, m = m, ev = ev, Xt = Xt, D = D, p = p, mruns = mruns)
            list(a = a, fit = fit, mom = mom, score = sum(diag(mom$Xi)))
        }, error = function(e) NULL)
        if (!is.null(obj) && is.finite(obj$score) && obj$score < score) {
            best <- obj
            score <- obj$score
        }
    }
    if (is.null(best))
        stop("No candidate has stable positive definite variance estimates")
    best
}

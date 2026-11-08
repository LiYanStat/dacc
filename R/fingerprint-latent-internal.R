# Internal calculations for the latent-factor framework of Li and Li (2026).
# W is a precision weight here; it is not the covariance weight used by Poly.

# Algorithm 1: two-step correction of control spikes and noise.
# Unlike inference, trial noise estimation for selecting K must not reject
# candidates solely because the last trial factor is at the BBP boundary.
.lf_adjusted <- function(values, k, n, p = length(values)) {
    gamma <- p / n
    noise.initial <- mean(values[seq.int(k + 1L, p)])

    if (!is.finite(noise.initial) || noise.initial <= 0) {
        stop("The trailing control spectrum must have positive average variance")
    }
    if (k == 0L) {
        return(list(spikes = numeric(), sigma_sq = noise.initial))
    }

    leading <- values[seq_len(k)]
    first <- .lf_spike_roots(leading, noise.initial, gamma)
    if (any(first <= 0)) {
        stop("The trial factor count gives nonpositive corrected spike strengths")
    }

    bias <- k + sum(noise.initial / first)
    noise <- noise.initial * (1 + gamma * bias / (p - k))
    spikes <- pmax(.lf_spike_roots(leading, noise, gamma), noise * sqrt(gamma))

    list(spikes = spikes, sigma_sq = noise)
}

# Inferential quantities from the corrected covariance. Inverse alignment is
# undefined at the detectability boundary even though noise remains estimable.
.lf_spikes <- function(values, n.factors, n.controls) {
    corrected <- .lf_adjusted(values, n.factors, n.controls)
    noise <- corrected$sigma_sq
    spikes <- corrected$spikes
    gamma <- length(values) / n.controls
    ratio <- spikes / noise
    alignment <- (1 - gamma / ratio^2) / (1 + gamma / ratio)

    if (any(!is.finite(alignment) | alignment <= sqrt(.Machine$double.eps))) {
        stop(paste(
            "Selected latent factors are at or below the detectability boundary;",
            "reduce n.factors. The variance formula requires positive alignment."
        ))
    }

    list(spikes = spikes, noise = noise, alignment = alignment)
}

# User-supplied PLY criterion and sequential KN test, evaluated on the full
# ordered spectrum of crossprod(controls) / n (including zeros when n < p).
.lf_select_K <- function(
    values,
    n,
    p,
    method = c("KN", "PLY"),
    conf_level = 0.95,
    delta = 0.05,
    K_max = min(n - 1L, floor(p * 0.1))
) {
    method <- match.arg(method)
    if (!.fp_integer_scalar(n, 2) || !.fp_integer_scalar(p, 1)) {
        stop("n and p must be positive integer control sample dimensions, n >= 2")
    }
    if (
        !is.numeric(values) || length(values) != p ||
        any(!is.finite(values) | values < 0) || any(diff(values) > 0)
    ) {
        stop("values must contain all p nonnegative control eigenvalues in decreasing order")
    }
    if (!.fp_integer_scalar(K_max, 0) || K_max >= min(n, p)) {
        stop("max.factors must be a nonnegative integer smaller than min(n, p)")
    }
    if (!.fp_positive(conf_level) || conf_level >= 1) {
        stop("factor.level must lie strictly between zero and one")
    }
    if (!.fp_positive(delta)) {
        stop("factor.delta must be positive")
    }
    if (K_max == 0L) {
        return(0L)
    }

    if (method == "PLY") {
        penalty <- (p / n + 2 * sqrt(p / n)) * (1 + n / p^(1 + delta)) / p
        counts <- seq.int(0L, K_max)
        noise <- vapply(counts, .lf_selection_noise, numeric(1L),
                        values = values, n = n, p = p)
        scores <- noise + counts * noise[K_max + 1L] * penalty
        return(as.integer(which.min(scores) - 1L))
    }

    if (!requireNamespace("RMTstat", quietly = TRUE)) {
        stop(paste(
            "Install RMTstat for KN spike selection,",
            "supply n.factors, or use factor.method = 'PLY'."
        ), call. = FALSE)
    }
    s <- RMTstat::qtw(conf_level, beta = 1)
    for (k in seq_len(K_max)) {
        d <- p - k
        mu <- (sqrt(n - 0.5) + sqrt(d - 0.5))^2 / n
        eta <- (sqrt(n - 0.5) + sqrt(d - 0.5)) *
            (1 / sqrt(n - 0.5) + 1 / sqrt(d - 0.5))^(1 / 3) / n
        noise <- .lf_selection_noise(k, values, n, p)
        if (values[k] <= noise * (mu + eta * s)) {
            return(as.integer(k - 1L))
        }
    }
    as.integer(K_max)
}

.lf_selection_noise <- function(k, values, n, p) {
    .lf_adjusted(values, k, n, p)$sigma_sq
}

# Larger real part of the quadratic roots, as specified in Algorithm 1.
.lf_spike_roots <- function(leading, noise, gamma) {
    b <- leading - (1 + gamma) * noise
    (b + sqrt(pmax(b^2 - 4 * gamma * noise^2, 0))) / 2
}

# Apply I + U diag(theta) U' without forming an N by N matrix.
.lf_weighted <- function(A, U, theta) {
    A + U %*% sweep(crossprod(U, A), 1L, theta, "*")
}

# Weighted TLS from the main-text criterion (7), including theta = -1.
.lf_tls <- function(theta, sigma.x2, Xtilde, Y, mruns, U) {
    N <- nrow(Xtilde)
    p <- ncol(Xtilde)
    D <- diag(1 / mruns, p)
    WX <- .lf_weighted(Xtilde, U, theta)
    WY <- .lf_weighted(matrix(Y, ncol = 1L), U, theta)
    G <- crossprod(Xtilde, WX) / N
    b <- crossprod(Xtilde, WY) / N
    A <- rbind(
        cbind(
            G * outer(sqrt(mruns), sqrt(mruns)),
            sqrt(sigma.x2 * mruns) * b
        ),
        c(
            sqrt(sigma.x2 * mruns) * b,
            sigma.x2 * sum(Y * WY) / N
        )
    )
    lambda <- min(eigen(A, symmetric = TRUE, only.values = TRUE)$values)
    beta <- as.vector(solve(G - lambda * D, b))
    residual <- Y - as.vector(Xtilde %*% beta)
    Wr <- .lf_weighted(matrix(residual, ncol = 1L), U, theta)

    list(
        beta = beta,
        residuals = residual,
        quadratic = as.numeric(crossprod(residual, Wr)) / N,
        scale = 1 + sigma.x2 * sum(beta^2 / mruns),
        gram = G
    )
}

# Trace approximations (10)--(13), using squared eigenvector alignments.
.lf_traces <- function(theta, spikes, noise, alignment, N) {
    K <- length(theta)
    delta <- noise + (sum(theta * (spikes * alignment + noise)) +
        sum(spikes)) / N
    phi <- sum(
        theta^2 * (alignment * spikes + noise)^2 +
        2 * theta * (alignment * spikes^2 +
            2 * noise * alignment * spikes + noise^2) +
        (spikes + noise)^2
    ) / N + (1 - K / N) * noise^2
    psi <- noise^2 * (1 + sum(2 * theta + theta^2) / N) +
        sum(spikes^2 * (alignment * theta + 1)^2) / N +
        noise * sum(spikes * (
            alignment * theta^2 + 2 * theta + 1 / alignment
        )) / N +
        noise * sum(spikes * (1 + alignment * (2 * theta + theta^2))) / N

    list(delta = delta, phi = phi, psi = psi)
}

# Equations (16)--(17), parameterized by variance s_x^2 for API consistency.
# An explicit inflation.grid is searched exactly; otherwise bracket the residual
# matching root and solve it numerically. No arbitrary fixed upper cap is hidden.
.lf_inflation <- function(theta, alpha, delta, Xtilde, Y, mruns, U, inflation.grid) {
    target <- alpha * delta
    discrepancy <- function(sigma.x2) {
        fit <- .lf_tls(theta, sigma.x2, Xtilde, Y, mruns, U)
        sigma.x2 * fit$quadratic / fit$scale - target
    }

    if (!is.null(inflation.grid)) {
        losses <- vapply(inflation.grid, function(s) {
            tryCatch(abs(discrepancy(s)), error = function(e) Inf)
        }, numeric(1L))
        i <- which.min(losses)
        if (!is.finite(losses[i]) || inflation.grid[i] <= 0) {
            stop("No positive finite inflation estimate on inflation.grid")
        }
        return(list(
            x = inflation.grid[i], z = inflation.grid[i] / alpha,
            boundary = i == length(inflation.grid)
        ))
    }

    upper <- 1
    value <- discrepancy(upper)
    while (value < 0 && upper < 2^20) {
        upper <- upper * 2
        value <- discrepancy(upper)
    }
    if (!is.finite(value) || value < 0) {
        stop("Could not bracket the inflation estimate; supply inflation.grid")
    }
    estimate <- stats::uniroot(
        discrepancy, interval = c(0, upper), tol = 1e-8
    )$root

    if (estimate <= 0) {
        stop("Inflation estimate must be positive")
    }
    list(x = estimate, z = estimate / alpha, boundary = FALSE)
}

# Re-estimate inflation for each candidate precision weight before scoring it.
.lf_evaluate <- function(theta, design) {
    trace <- .lf_traces(
        theta, design$spikes, design$noise, design$alignment, design$N
    )
    if (trace$delta <= 0 || trace$phi <= 0) {
        stop("Nonpositive latent trace estimates")
    }
    inflation <- if (design$estimate.inflation) {
        .lf_inflation(
            theta, design$alpha, trace$delta, design$Xtilde, design$Y,
            design$mruns, design$U, design$inflation.grid
        )
    } else {
        list(x = 1, z = 1, boundary = FALSE)
    }
    fit <- .lf_tls(
        theta, inflation$x, design$Xtilde, design$Y, design$mruns, design$U
    )
    variance <- .lf_variance(theta, inflation, trace, fit, design)

    list(
        theta = theta, inflation = inflation, trace = trace,
        fit = fit, Xi = variance, objective = sum(diag(variance))
    )
}

# Invalid finite-sample candidates are excluded from the variance minimization.
.lf_objective <- function(theta, design) {
    tryCatch(
        .lf_evaluate(theta, design)$objective,
        error = function(e) 1e100
    )
}

# Local constrained optimization from naive precision and identity starts.
.lf_select <- function(design, theta, optim.control) {
    if (!is.null(theta)) {
        return(.lf_evaluate(theta, design))
    }
    K <- length(design$spikes)
    if (K == 0L) {
        return(.lf_evaluate(numeric(), design))
    }
    starts <- list(
        -design$spikes / (design$spikes + design$noise),
        rep(0, K)
    )
    best <- NULL
    converged <- FALSE
    for (start in starts) {
        # Preserve feasible starts even if the optimizer encounters an invalid
        # variance estimate at a later candidate.
        candidate <- tryCatch(
            .lf_evaluate(start, design), error = function(e) NULL
        )
        if (!is.null(candidate) &&
            (is.null(best) || candidate$objective < best$objective)) {
            best <- candidate
        }
        result <- tryCatch(
            stats::optim(
                start, .lf_objective, design = design,
                method = "L-BFGS-B", lower = rep(-1, K),
                control = optim.control
            ),
            error = function(e) NULL
        )
        if (is.null(result)) next
        candidate <- tryCatch(
            .lf_evaluate(result$par, design), error = function(e) NULL
        )
        if (!is.null(candidate) && result$convergence == 0) {
            converged <- TRUE
        }
        if (!is.null(candidate) &&
            (is.null(best) || candidate$objective <= best$objective)) {
            candidate$convergence <- result$convergence
            best <- candidate
        }
    }
    if (is.null(best)) {
        stop("No candidate has a stable positive definite latent variance estimate")
    }
    if (!converged) {
        warning("Latent weight optimization did not converge; returning the best feasible candidate")
    }
    best
}

# Unbiased within-forcing variation divided by average control variation.
.lf_alpha <- function(Xtilde, mruns, ensembles, control.variation) {
    N <- nrow(Xtilde)
    p <- ncol(Xtilde)
    if (!is.list(ensembles) || length(ensembles) != p || sum(mruns) <= p) {
        stop("Factor estimation requires individual ensembles and within-ensemble degrees of freedom")
    }
    ss <- 0
    for (i in seq_len(p)) {
        E <- ensembles[[i]]
        if (!.fp_finite_matrix(E) ||
            !identical(dim(E), c(as.integer(mruns[i]), N))) {
            stop("Each ensemble must be an mruns[i] by N finite numeric matrix")
        }
        if (!isTRUE(all.equal(
            as.numeric(colMeans(E)), as.numeric(Xtilde[, i]), tolerance = 1e-7
        ))) {
            stop("Ensemble means disagree with Xtilde")
        }
        ss <- ss + sum(sweep(E, 2L, Xtilde[, i])^2)
    }
    ratio <- ss / (N * (sum(mruns) - p)) / control.variation
    if (!is.finite(ratio) || ratio <= 0) {
        stop("Positive within-ensemble and control variation are required")
    }
    ratio
}

# Theorem 4, with G and Omega from equations (14)--(15).
.lf_variance <- function(theta, inflation, trace, fit, design) {
    X <- design$Xtilde
    p <- ncol(X)
    D <- diag(1 / design$mruns, p)
    UTX <- crossprod(design$U, X)
    W2X <- .lf_weighted(X, design$U, 2 * theta + theta^2)
    spike.coefficient <- design$spikes / design$alignment *
        (design$alignment * theta + 1)^2
    WtildeX <- design$U %*% sweep(UTX, 1L, spike.coefficient, "*") +
        design$noise * W2X

    G <- fit$gram - design$alpha * trace$delta * D
    Omega <- crossprod(X, WtildeX) / design$N -
        design$alpha * trace$psi * D
    G <- (G + t(G)) / 2
    Omega <- (Omega + t(Omega)) / 2
    if (any(!is.finite(G)) || any(!is.finite(Omega))) {
        stop("Nonfinite latent covariance estimates")
    }
    chol(G)
    chol(Omega)
    inverse <- solve(G)
    correction <- solve(
        diag(design$mruns, p) + inflation$x * tcrossprod(fit$beta)
    )
    Xi <- fit$scale / inflation$z * inverse %*% (
        Omega + design$alpha * trace$phi * correction
    ) %*% inverse
    Xi <- (Xi + t(Xi)) / 2
    chol(Xi)
    Xi
}

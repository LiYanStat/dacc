if (file.exists(file.path("R", "fingerprintLatent.R"))) {
    for (file in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(file)
} else {
    library(dacc)
    for (name in c(".lf_spikes", ".lf_spike_roots", ".lf_traces")) {
        assign(name, getFromNamespace(name, "dacc"))
    }
}

# Trace formulas (10)--(13) against direct matrices with known alignments.
set.seed(15)
N <- 30L
basis <- qr.Q(qr(matrix(rnorm(N * 4), N)))
true.U <- basis[, 1:2]
z <- c(0.8, 0.6)
empirical.U <- sweep(true.U, 2, sqrt(z), "*") +
    sweep(basis[, 3:4], 2, sqrt(1 - z), "*")
pi <- c(7, 4)
noise <- 1.3
theta <- c(-1, 0.4)
C <- noise * diag(N) + true.U %*% diag(pi) %*% t(true.U)
W <- diag(N) + empirical.U %*% diag(theta) %*% t(empirical.U)
Wt <- empirical.U %*% diag(pi / z * (z * theta + 1)^2) %*%
    t(empirical.U) + noise * W %*% W
tr <- .lf_traces(theta, pi, noise, z, N)
stopifnot(
    abs(tr$delta - sum(diag(W %*% C)) / N) < 1e-12,
    abs(tr$phi - sum(diag(W %*% C %*% W %*% C)) / N) < 1e-12,
    abs(tr$psi - sum(diag(Wt %*% C)) / N) < 1e-12
)

# Algorithm 1, including the noise bias adjustment.
known.pi <- c(15, 8)
gamma <- 0.5
noise0 <- 2
leading <- noise0 + known.pi + gamma * noise0 * (known.pi + noise0) / known.pi
stopifnot(max(abs(.lf_spike_roots(leading, noise0, gamma) - known.pi)) < 1e-12)
corrected <- .lf_spikes(c(leading, rep(noise0, 38)), 2, 80)
expected.noise <- noise0 * (1 + gamma * (2 + sum(noise0 / known.pi)) / 38)
stopifnot(abs(corrected$noise - expected.noise) < 1e-12)

# Detectable factors, unequal ensemble sizes, and variance heterogeneity.
set.seed(95)
N <- 120L
m <- 180L
K <- 2L
sizes <- c(25L, 35L)
U <- qr.Q(qr(matrix(rnorm(N * K), N)))
Sigma <- diag(N) + U %*% diag(c(18, 12)) %*% t(U)
L <- chol(Sigma)
truth <- cbind(3 * sin(seq_len(N) / 8), 4 * cos(seq_len(N) / 11))
Z <- matrix(rnorm(m * N), m) %*% L * sqrt(1.5)
ensembles <- lapply(seq_len(2), function(i) {
    sweep(matrix(rnorm(sizes[i] * N), sizes[i]) %*% L * sqrt(1.8),
          2, truth[, i], "+")
})
X <- vapply(ensembles, colMeans, numeric(N))
Y <- as.vector(truth %*% c(1, 0.7)) + as.vector(rnorm(N) %*% L)
fixed <- fingerprintLatent(
    X, Y, sizes, Z, n.factors = K,
    estimate.inflation = FALSE, theta = c(-0.7, -1)
)
stopifnot(
    inherits(fixed, "fingerprint"), inherits(fixed, "fingerprintLatent"),
    fixed$sigma.x2 == 1, fixed$sigma.z2 == 1,
    !any(c("Q", "V", "history", "settings", "diagnostic") %in% names(fixed)),
    all(is.finite(fixed$ci))
)

# Independent precision-weighted SVD solution, including a removed direction.
e <- eigen(fixed$weight, symmetric = TRUE)
A <- sweep(crossprod(e$vectors, cbind(sweep(X, 2, sqrt(sizes), "*"), Y)),
           1, sqrt(pmax(e$values, 0)), "*")
v <- svd(A)$v[, 3]
reference <- -sqrt(sizes) * v[1:2] / v[3]
stopifnot(max(abs(reference - fixed$beta)) < 1e-8)

# Theorem 4 against a direct matrix calculation, including D in Omega.
S <- crossprod(Z) / m
e <- eigen(S, symmetric = TRUE)
sp <- .lf_spikes(pmax(e$values, 0), K, m)
Ut <- e$vectors[, seq_len(K), drop = FALSE]
t <- fixed$theta
tr <- .lf_traces(t, sp$spikes, sp$noise, sp$alignment, N)
D <- diag(1 / sizes)
G <- crossprod(X, fixed$weight %*% X) / N - tr$delta * D
Wt <- Ut %*% diag(sp$spikes / sp$alignment * (sp$alignment * t + 1)^2) %*%
    t(Ut) + sp$noise * fixed$weight %*% fixed$weight
Omega <- crossprod(X, Wt %*% X) / N - tr$psi * D
scale <- 1 + sum(fixed$beta^2 / sizes)
Xi <- scale * solve(G) %*% (
    Omega + tr$phi * solve(diag(sizes) + tcrossprod(fixed$beta))
) %*% solve(G)
stopifnot(max(abs(fixed$var - Xi / N)) < 1e-10)
Q <- as.numeric(crossprod(fixed$residuals, fixed$weight %*% fixed$residuals)) / N -
    scale * tr$delta
V <- 2 / N * scale^2 * tr$phi
stopifnot(
    abs(fixed$diagnostic.stat - Q / sqrt(V)) < 1e-10,
    max(abs(fixed$standardized.residuals - fixed$residuals / sqrt(diag(S)))) < 1e-12
)

# Factor estimation and weight selection; exact matching degenerates the test.
fit <- fingerprintLatent(X, Y, sizes, Z, K, ensembles = ensembles)
within <- sum(vapply(seq_len(2), function(i) {
    sum(sweep(ensembles[[i]], 2, X[, i])^2)
}, numeric(1))) / (N * (sum(sizes) - 2))
stopifnot(
    fit$sigma.x2 > 0, fit$sigma.z2 > 0,
    all(is.finite(fit$ci)), all(fit$theta >= -1),
    abs(fit$p.value - 0.5) < 1e-6,
    abs(fit$sigma.x2 / fit$sigma.z2 - within / mean(diag(S))) < 1e-10
)
dispatched <- fingerprint(X, Y, sizes, Z, method = "latent", n.factors = K,
                          ensembles = ensembles)
stopifnot(identical(fit$beta, dispatched$beta), identical(fit$var, dispatched$var))

# Explicit grids, the no-spike special case, and singular control covariance.
grid <- seq(0, 5, length.out = 101)
grid.fit <- fingerprintLatent(X, Y, sizes, Z, K, ensembles = ensembles,
                              theta = c(-0.7, -0.8), inflation.grid = grid)
stopifnot(grid.fit$sigma.x2 %in% grid)
isotropic <- fingerprintLatent(X, Y, sizes, Z, 0, estimate.inflation = FALSE)
stopifnot(length(isotropic$theta) == 0L, identical(isotropic$weight, diag(N)))
for (estimate in c(FALSE, TRUE)) {
    singular <- fingerprintLatent(X, Y, sizes, Z[1:80, ], K,
        ensembles = ensembles, estimate.inflation = estimate, theta = c(-0.7, -0.8))
    stopifnot(all(is.finite(singular$ci)))
}
stopifnot(
    inherits(try(fingerprintLatent(X, Y, sizes, Z), silent = TRUE), "try-error"),
    inherits(try(fingerprintLatent(X, Y, sizes, Z, K), silent = TRUE), "try-error"),
    inherits(try(fingerprintLatent(X, Y, sizes, Z, K,
        theta = c(-2, 0)), silent = TRUE), "try-error"),
    inherits(try(.lf_spikes(rep(1, 20), 1, 40), silent = TRUE), "try-error")
)
file <- tempfile(fileext = ".png")
grDevices::png(file)
visible <- withVisible(plot(fit))
grDevices::dev.off()
stopifnot(!visible$visible, identical(visible$value, fit), file.info(file)$size > 0)
cat("Latent framework checks passed\n")

if (file.exists(file.path("R", "fingerprint.R"))) {
    for (file in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
        source(file)
    }
} else {
    library(dacc)
}

# Test routing independently of the numerical implementations.
routing <- new.env(parent = environment(fingerprint))
routing$fingerprintPoly <- function(...) list(backend = "poly", arguments = list(...))
routing$fingerprintLatent <- function(...) list(backend = "latent", arguments = list(...))
routing$fingerprintLegacy <- function(...) list(backend = "legacy", arguments = list(...))
route <- fingerprint
environment(route) <- routing
X <- matrix(seq_len(20), 10, 2)
Y <- seq_len(10)
Z <- matrix(0, 30, 10)

fit <- route(X, Y, c(3, 4), Z, degree = 3, estimate.inflation = FALSE)
stopifnot(fit$backend == "poly", fit$arguments$degree == 3)
fit <- route(X, Y, c(3, 4), Z, method = "latent", n.factors = 2)
stopifnot(fit$backend == "latent", fit$arguments$n.factors == 2)
for (method in c("EE", "PBC", "TS", "TS.TempSt")) {
    fit <- route(X, Y, c(3, 4), Z, method = method, S = 2, T = 5)
    stopifnot(
        fit$backend == "legacy", fit$arguments$method == method,
        identical(fit$arguments$ctlruns.sigma, Z),
        identical(fit$arguments$ctlruns.bhvar, Z)
    )
}
stopifnot(
    identical(eval(formals(fingerprint)$method),
        c("poly", "linear", "latent", "EE", "PBC", "TS", "TS.TempSt")),
    !"legacy.method" %in% names(formals(fingerprint)),
    inherits(try(route(X, Y, c(3, 4), Z, method = "legacy"),
        silent = TRUE), "try-error")
)

# The explicit linear method fixes degree one without fixing inflation coefficients.
for (estimate in c(TRUE, FALSE)) {
    fit <- route(X, Y, c(3, 4), Z, method = "linear",
                 estimate.inflation = estimate)
    stopifnot(fit$backend == "poly", fit$arguments$degree == 1L,
              identical(fit$arguments$estimate.inflation, estimate))
}
for (degree in list(3, 0, NA_real_, NULL, "1", c(1, 3))) {
    fit <- route(X, Y, c(3, 4), Z, method = "linear", degree = degree)
    stopifnot(fit$backend == "poly", identical(fit$arguments$degree, 1L))
}
fit <- route(X, Y, c(3, 4), Z, method = "linear", degree = 1)
stopifnot(fit$arguments$degree == 1L)

# Renaming the polynomial backend must not change its estimates.
set.seed(89)
X <- cbind(3 * sin(seq_len(40) / 4), 3 * cos(seq_len(40) / 5))
Y <- as.vector(X %*% c(1, 0.5) + rnorm(40))
Z <- matrix(rnorm(120 * 40), 120)
a <- fingerprint(X, Y, c(30, 30), Z, degree = 1, estimate.inflation = FALSE)
b <- fingerprintPoly(X, Y, c(30, 30), Z, degree = 1, estimate.inflation = FALSE)
stopifnot(
    inherits(a, "fingerprintPoly"), inherits(a, "fingerprint"),
    identical(a$beta, b$beta), identical(a$var, b$var), identical(a$ci, b$ci)
)

linear <- fingerprint(X, Y, c(30, 30), Z,
                      method = "linear", degree = 3, estimate.inflation = FALSE)
stopifnot(identical(linear$beta, b$beta), identical(linear$var, b$var),
          identical(linear$ci, b$ci))

# Original estimator and bootstrap intervals remain unchanged by dispatch.
Y <- matrix(Y, ncol = 1)
set.seed(55)
a <- fingerprint(
    X, Y, c(30, 30), Z, method = "TS",
    S = 4, T = 10, B = 5, cov.method = "l2"
)
set.seed(55)
b <- fingerprintLegacy(
    X, Y, c(30, 30), Z, Z, 4, 10,
    B = 5, method = "TS", cov.method = "l2"
)
stopifnot(identical(a, b))
# A separate interval control sample must also pass through unchanged.
Z2 <- Z + 1
fit <- route(X, Y, c(30, 30), Z, method = "PBC", ctlruns.bhvar = Z2)
stopifnot(identical(fit$arguments$ctlruns.bhvar, Z2))
cat("Common entry point and legacy equivalence checks passed\n")

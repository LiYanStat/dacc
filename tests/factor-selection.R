if (file.exists(file.path("R", "fingerprintLatent.R"))) {
    for (file in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(file)
} else {
    library(dacc)
    for (name in c(".lf_adjusted", ".lf_select_K")) {
        assign(name, getFromNamespace(name, "dacc"))
    }
}

# Original user-supplied selector, with the manuscript noise correction supplied.
reference <- function(values, n, p, method, conf_level = 0.95,
                      delta = 0.05, K_max = min(n - 1, floor(p * 0.1))) {
    if (K_max == 0L) return(0L)
    noise <- function(k) .lf_adjusted(values, k, n, p)$sigma_sq
    if (method == "PLY") {
        penalty <- (p / n + 2 * sqrt(p / n)) * (1 + n / p^(1 + delta)) / p
        scores <- vapply(0:K_max, noise, numeric(1)) +
            (0:K_max) * noise(K_max) * penalty
        return(which.min(scores) - 1L)
    }
    s <- RMTstat::qtw(conf_level, beta = 1)
    for (k in seq_len(K_max)) {
        d <- p - k
        mu <- (sqrt(n - 0.5) + sqrt(d - 0.5))^2 / n
        eta <- (sqrt(n - 0.5) + sqrt(d - 0.5)) *
            (1 / sqrt(n - 0.5) + 1 / sqrt(d - 0.5))^(1 / 3) / n
        if (values[k] <= noise(k) * (mu + eta * s)) return(k - 1L)
    }
    K_max
}

set.seed(95)
N <- 120L
m <- 180L
U <- qr.Q(qr(matrix(rnorm(N * 2), N)))
L <- chol(diag(N) + U %*% diag(c(18, 12)) %*% t(U))
Z <- matrix(rnorm(m * N), m) %*% L
values <- pmax(eigen(crossprod(Z) / m, symmetric = TRUE)$values, 0)
methods <- "PLY"
if (requireNamespace("RMTstat", quietly = TRUE)) methods <- c(methods, "KN")
for (method in methods) {
    for (controls in list(Z, Z[1:80, ])) {
        n <- nrow(controls)
        spectrum <- pmax(eigen(crossprod(controls) / n, symmetric = TRUE)$values, 0)
        k <- .lf_select_K(spectrum, n, N, method)
        stopifnot(
            identical(k, as.integer(reference(spectrum, n, N, method))),
            identical(k, .lf_select_K(37 * spectrum, n, N, method)),
            identical(.lf_select_K(spectrum, n, N, method, K_max = 1),
                      as.integer(reference(spectrum, n, N, method, K_max = 1)))
        )
    }
    stopifnot(.lf_select_K(values, m, N, method) == 2L)
}

# The supplied default cap gives no candidate factors when dimension < 10.
stopifnot(identical(.lf_select_K(rep(1, 8), 20, 8, "KN"), 0L))
stopifnot(identical(.lf_select_K(values, m, N, "PLY", K_max = 0), 0L))

# No-spike controls should allow zero factors, rather than force at least one.
set.seed(91)
noise <- matrix(rnorm(m * N), m)
noise.values <- pmax(eigen(crossprod(noise) / m, symmetric = TRUE)$values, 0)
for (method in methods) {
    stopifnot(.lf_select_K(noise.values, m, N, method) == 0L)
}

# Missing optional dependency has an actionable error; PLY remains usable.
environment.copy <- new.env(parent = environment(.lf_select_K))
environment.copy$requireNamespace <- function(...) FALSE
without.RMTstat <- .lf_select_K
environment(without.RMTstat) <- environment.copy
error <- tryCatch(without.RMTstat(values, m, N, "KN"), error = conditionMessage)
stopifnot(is.character(error), grepl("RMTstat", error),
          without.RMTstat(values, m, N, "PLY") == 2L)

# Auto-selection and an explicit equal count must give identical fit numerics.
X <- cbind(3 * sin(seq_len(N) / 8), 4 * cos(seq_len(N) / 11))
Y <- as.vector(X %*% c(1, 0.7) + rnorm(N))
for (method in methods) {
    automatic <- fingerprint(X, Y, c(25, 35), Z, method = "latent",
                             factor.method = method, estimate.inflation = FALSE)
    explicit <- fingerprintLatent(X, Y, c(25, 35), Z, n.factors = 2,
                                  estimate.inflation = FALSE)
    stopifnot(
        automatic$n.factors == 2L,
        identical(automatic$beta, explicit$beta),
        identical(automatic$var, explicit$var),
        identical(automatic$theta, explicit$theta)
    )
}
override <- fingerprintLatent(X, Y, c(25, 35), Z,
                              n.factors = 0, estimate.inflation = FALSE)
stopifnot(override$n.factors == 0L, length(override$theta) == 0L)

# Validate selector inputs before indexing spectra or taking square roots.
for (bad in list(-1, N, 1.5, NA_real_)) {
    stopifnot(inherits(try(.lf_select_K(values, m, N, K_max = bad),
                          silent = TRUE), "try-error"))
}
stopifnot(
    inherits(try(.lf_select_K(values[-1], m, N), silent = TRUE), "try-error"),
    inherits(try(.lf_select_K(rev(values), m, N), silent = TRUE), "try-error"),
    inherits(try(.lf_select_K(values, m, N, conf_level = 1), silent = TRUE), "try-error")
)
cat("Automatic PLY and KN factor-selection checks passed\n")

# KN is the default in both the fitting interface and the internal selector.
stopifnot(
    identical(eval(formals(fingerprintLatent)$factor.method), c("KN", "PLY")),
    identical(eval(formals(.lf_select_K)$method), c("KN", "PLY"))
)
if (requireNamespace("RMTstat", quietly = TRUE)) {
    stopifnot(identical(.lf_select_K(values, m, N),
                        .lf_select_K(values, m, N, method = "KN")))
    default <- fingerprintLatent(X, Y, c(25, 35), Z,
                                 estimate.inflation = FALSE)
    explicit.KN <- fingerprintLatent(X, Y, c(25, 35), Z,
                                     estimate.inflation = FALSE,
                                     factor.method = "KN")
    default$call <- explicit.KN$call <- NULL
    stopifnot(identical(default, explicit.KN))
}

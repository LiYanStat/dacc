if (file.exists(file.path("R", "fingerprint.R"))) {
    for (file in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(file)
} else {
    library(dacc)
}

set.seed(2026)
N <- 80
n <- 30
m <- 160
u <- rep(1 / sqrt(N), N)
signal <- 3 * sin(seq_len(N) / 6)
Z <- matrix(rnorm(m * N), m) + outer(rnorm(m), sqrt(20) * u)
E <- matrix(rnorm(n * N), n) + outer(rnorm(n), sqrt(20) * u)
E <- sweep(E, 2, signal, "+")
X <- matrix(colMeans(E), ncol = 1)
Y <- signal + rnorm(N) + sqrt(20) * rnorm(1) * u

for (method in c("poly", "linear", "latent")) {
    specific <- if (method == "poly") list(degree = 1, grid.size = 5) else
        if (method == "latent") list(n.factors = 1, theta = -0.7) else list()
    common <- c(list(Xtilde = X, Y = Y, mruns = n, ctlruns = Z,
                     ensembles = list(E), method = method), specific)
    for (estimate in c(FALSE, TRUE)) {
        canonical <- do.call(fingerprint, c(common,
            list(estimate.inflation = estimate,
                 inflation.grid = seq(0, 8, length.out = 81))))
        stopifnot(all(is.finite(canonical$ci)))
        if (!estimate) {
            stopifnot(canonical$sigma.x2 == 1, canonical$sigma.z2 == 1)
        }
    }
    for (removed in list(list(estimate.factors = FALSE), list(factor.grid = 0:3))) {
        error <- tryCatch(do.call(fingerprint, c(common, removed)),
                          error = conditionMessage)
        stopifnot(is.character(error), grepl("unused argument", error))
    }
}
# Latent factor-count controls keep their own names and defaults.
stopifnot(
    "n.factors" %in% names(formals(fingerprintLatent)),
    "factor.method" %in% names(formals(fingerprintLatent)),
    "estimate.inflation" %in% names(formals(fingerprintLatent)),
    "inflation.grid" %in% names(formals(fingerprintPoly))
)
for (fun in list(fingerprintPoly, fingerprintLatent)) {
    stopifnot(!any(c("estimate.factors", "factor.grid") %in% names(formals(fun))))
}
cat("Variance inflation argument checks passed\n")

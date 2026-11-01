if (file.exists(file.path("R", "fingerprint.R"))) {
  source(file.path("R", "fingerprint-internal.R"))
  source(file.path("R", "fingerprint.R"))
} else library(dacc)
set.seed(403)
N <- 100; sizes <- c(30L, 35L)
X <- cbind(sin(seq_len(N)/7)*3, cos(seq_len(N)/9)*3)
E <- lapply(seq_len(2), function(i) sweep(matrix(rnorm(sizes[i]*N), sizes[i]), 2, X[, i], "+"))
Xm <- vapply(E, colMeans, numeric(N))
Z <- matrix(rnorm(200*N), 200)
y <- as.vector(X %*% c(1, 0.7) + rnorm(N))
fixed <- fingerprint(Xm, y, sizes, Z, degree = 1, estimate.factors = FALSE)
stopifnot(fixed$sigma.x2 == 1, fixed$sigma.z2 == 1, !any(c("history", "settings", "diagnostic") %in% names(fixed)),
          length(fixed$coefficients) == 2, all(is.finite(fixed$ci)))
# Independent generalized TLS solution via a whitened augmented SVD.
W <- eigen(fixed$weight, symmetric = TRUE)
A <- crossprod(W$vectors, cbind(sweep(Xm, 2, sqrt(sizes), "*"), y))/sqrt(W$values)
v <- svd(A)$v[, 3]
reference <- -sqrt(sizes)*v[1:2]/v[3]
stopifnot(max(abs(reference-fixed$beta)) < 1e-8)
# Spectral trace/derivative identities for the degree-one partial fractions.
a <- fixed$coefficients[2]; S <- crossprod(Z)/nrow(Z); r <- -1/a
h <- 1/(eigen(S, symmetric = TRUE)$values-r); rho <- N/nrow(Z)
th <- 1/(1-rho-rho*r*mean(h))
pi <- (th-1)/rho
pip <- (mean(h)+r*mean(h^2))*th^2
K <- th^2*(r*pip+pi)/a^2
scale <- 1+sum(fixed$beta^2/sizes)
reference.V <- 2/N*scale^2*K
reference.Q <- as.numeric(crossprod(fixed$residuals,
  solve(fixed$weight, fixed$residuals)))/N - scale*pi/a
stopifnot(abs(fixed$diagnostic.stat - reference.Q/sqrt(reference.V)) < 1e-10)
for (k in c(1L, 2L, 3L)) {
  fit <- fingerprint(Xm, y, sizes, Z, E, degree = k, grid.size = 5,
                        factor.grid = seq(0, 4, length.out = 101))
  stopifnot(length(fit$coefficients) == k+1, !any(c("history", "settings", "diagnostic") %in% names(fit)),
            fit$sigma.x2 > 0, fit$sigma.z2 > 0,
            all(eigen(fit$var, symmetric = TRUE)$values > 0),
            fit$p.value >= 0, fit$p.value <= 1)
}
# Fixed factors require no ensembles, including with default cubic degree.
cubic <- fingerprint(Xm, y, sizes, Z, estimate.factors = FALSE)
stopifnot(length(cubic$coefficients) == 4, cubic$sigma.x2 == 1)
# Missing ensembles and an inadmissible decreasing polynomial are rejected.
stopifnot(inherits(try(fingerprint(Xm, y, sizes, Z), silent = TRUE), "try-error"),
          inherits(try(fingerprint(Xm, y, sizes, Z, degree = 1,
            estimate.factors = FALSE, coefficient.grid = matrix(-1)), silent = TRUE), "try-error"))
cat("Polynomial fingerprinting tests passed\n")
# Large TLS bounds must not hide order-one factors on a coarse grid.
set.seed(12)
n <- 60L
truth <- 3*sin(seq_len(n)/6)
controls <- matrix(rnorm(120*n), 120)
observed <- truth+rnorm(n)
members <- list(sweep(matrix(rnorm(30*n),30), 2, truth, "+"))
means <- matrix(colMeans(members[[1]]), ncol = 1)
adaptive <- fingerprint(means, observed, 30, controls, members, grid.size = 5)
stopifnot(adaptive$sigma.x2 > 0, adaptive$sigma.x2 < 10)
singular <- fingerprint(means, observed, 30, controls[1:40, ],
                           estimate.factors = FALSE, grid.size = 5)
stopifnot(all(is.finite(singular$ci)))

stopifnot("diagnostic.stat" %in% names(fixed),
          !any(c("statistic", "Q", "V") %in% names(fixed)),
          isTRUE(all.equal(fixed$p.value,
            pnorm(fixed$diagnostic.stat, lower.tail = FALSE))))

# Both inference branches must support singular S when m < N.
set.seed(403)
N <- 100L; sizes <- c(30L, 35L)
X <- cbind(3*sin(seq_len(N)/7), 3*cos(seq_len(N)/9))
E <- lapply(seq_len(2), function(i)
  sweep(matrix(rnorm(sizes[i]*N), sizes[i]), 2, X[, i], "+"))
Xm <- vapply(E, colMeans, numeric(N))
Z <- matrix(rnorm(40*N), 40)
y <- as.vector(X %*% c(1, 0.7) + rnorm(N))
stopifnot(nrow(Z) < N, qr(crossprod(Z))$rank < N)
for (k in c(1L, 3L)) for (estimate in c(FALSE, TRUE)) {
  fit <- fingerprint(Xm, y, sizes, Z, ensembles = E, degree = k,
                     estimate.factors = estimate, grid.size = 5)
  stopifnot(all(is.finite(fit$ci)), is.finite(fit$p.value),
            min(eigen(fit$weight, symmetric = TRUE)$values) > 0)
}
cat("Singular control covariance regression checks passed\n")

# Exercise S3 plot dispatch and custom coordinates on an actual graphics device.
if (file.exists(file.path("R", "plot.fingerprint.R"))) {
  source(file.path("R", "plot.fingerprint.R"))
}
stopifnot(inherits(fixed, "fingerprint"))
plot.file <- tempfile(fileext = ".png")
grDevices::png(plot.file, width = 900, height = 600)
result <- withVisible(plot(fixed))
stopifnot(!result$visible, identical(result$value, fixed))
plot(fixed, index = seq_along(fixed$residuals) + 1900,
     xlab = "Year", col = "steelblue", type = "b")
stopifnot(inherits(try(plot(fixed, index = 1:3), silent = TRUE), "try-error"))
grDevices::dev.off()
stopifnot(file.info(plot.file)$size > 0)
cat("Residual plot checks passed\n")

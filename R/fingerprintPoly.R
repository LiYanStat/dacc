##' Detection and attribution with polynomial shrinkage under variance heterogeneity
##'
##' Implements the method proposed in work by Li and Li (2026) for detection
##' and attribution (D&A) analysis of climate change using polynomial shrinkage
##' under variance heterogeneity. The framework allows forced simulations and
##' control runs to have different internal variability magnitudes from observations,
##' while sharing a common covariance structure. It estimates forcing scaling factors,
##' optionally estimates variance inflation coefficients, and provides asymptotic confidence
##' intervals and approximate residual diagnostics.
##'
##' @details
##' \subsection{Model}{
##' Li and Li (2026) extend optimal fingerprinting to allow different magnitudes
##' of internal variability in observations, forced simulations and control runs:
##' \deqn{Y=X\beta+\epsilon,\qquad \epsilon\sim N(0,\Sigma),}
##' \deqn{\widetilde X_{ik}=X_i+\eta_{ik},\qquad
##' \eta_{ik}\sim N(0,\sigma_x^2\Sigma),\quad k=1,\ldots,n_i,}
##' \deqn{Z_j\sim N(0,\sigma_z^2\Sigma),\quad j=1,\ldots,m.}
##' Here \eqn{Y\in\mathbb{R}^N} is the observed climate response,
##' \eqn{X=(X_1,\ldots,X_p)} contains the true responses to \eqn{p} forcings,
##' and \eqn{\beta} contains their scaling factors. The errors are mutually
##' independent and share the covariance structure \eqn{\Sigma}. The coefficients
##' \eqn{\sigma_x^2} and \eqn{\sigma_z^2} measure simulation variability relative
##' to observations; setting both to one gives the conventional equal-variability model.
##'
##' The inputs \code{Xtilde} and \code{mruns} specify the ensemble means
##' \eqn{\widetilde X_i=n_i^{-1}\sum_k\widetilde X_{ik}} and ensemble sizes
##' \eqn{n_i}. Write \eqn{D=\operatorname{diag}(1/n_1,\ldots,1/n_p)}.
##' }
##' \subsection{Polynomial shrinkage and scaling factors}{
##' The control covariance \eqn{S=m^{-1}\sum_j Z_jZ_j^\top} may be singular,
##' including when \eqn{m<N}. It is regularized using
##' \deqn{W=f(S),\qquad f(t)=1+\sum_{j=1}^k\ell_jt^j.}
##' Since \eqn{f(0)=1} and \eqn{f} is nondecreasing on the sample spectrum,
##' \eqn{W} is positive definite even when \eqn{S} is singular.
##' The degree \eqn{k} is controlled by \code{degree}: one gives linear shrinkage,
##' while three is recommended in Li and Li (2026). Candidate polynomials must be
##' nondecreasing on \eqn{[0,\psi_1]}, satisfy \eqn{f(\psi_1)\leq\psi_2}, and
##' have roots separated by at least \eqn{\psi_3}.
##'
##' For a given \eqn{W} and \eqn{s=\sigma_x^2}, the scaling factors minimize
##' the weighted total least squares criterion
##' \deqn{\widehat\beta(s,W)=\arg\min_\beta
##' \frac{(Y-\widetilde X\beta)^\top W^{-1}(Y-\widetilde X\beta)}
##' {1+s\beta^\top D\beta}.}
##' The polynomial coefficients are selected by minimizing the estimated total
##' asymptotic variance \eqn{\operatorname{tr}\widehat\Xi(f)} over an admissible
##' grid. The default uses 21 points per coefficient; \code{grid.size} and
##' \code{coefficient.grid} control the search resolution. This is a finite-grid
##' approximation to the optimization in Li and Li (2026).
##' }
##' \subsection{Estimation of variance inflation coefficients and inference}{
##' The method estimates both variance inflation coefficients, \eqn{\sigma_x^2} and
##' \eqn{\sigma_z^2}, to account for differences in internal variability between
##' observations, forced simulations and control runs. This is the default behavior
##' (\code{estimate.inflation = TRUE}). Estimation uses the individual forced runs
##' in \code{ensembles}, in addition to the ensemble means and control runs.
##'
##' The within-ensemble variation estimates the marginal variability of the
##' forced simulations:
##' \deqn{V_X=\frac{\sum_{i=1}^p\sum_{k=1}^{n_i}
##' \|\widetilde X_{ik}-\widetilde X_i\|^2}{N(\sum_i n_i-p)}.}
##' For a given weight matrix \eqn{W}, the forced-simulation inflation coefficient is
##' estimated by matching this variation to the residual variation implied by
##' candidate values of \eqn{\sigma_x^2}:
##' \deqn{\widehat\sigma_x^2=\arg\min_{s\in\mathcal G}
##' \left|\frac{s\|Y-\widetilde X\widehat\beta(s,W)\|^2}
##' {N\{1+s\widehat\beta(s,W)^\top D\widehat\beta(s,W)\}}-V_X\right|,}
##' where \eqn{\mathcal G} is the grid supplied through \code{inflation.grid} or
##' constructed automatically. The control-run inflation coefficient is then estimated
##' from the relative variation in controls and forced simulations:
##' \deqn{\widehat\sigma_z^2=\widehat\sigma_x^2\frac{V_Z}{V_X},
##' \qquad V_Z=\operatorname{tr}(S)/N.}
##'
##' Initial estimates of both coefficients are obtained using normalized Ledoit--Wolf
##' weights. These estimates are used to select the polynomial shrinkage rule.
##' The coefficients are then re-estimated using the selected weight matrix, and the
##' polynomial is reselected using the updated coefficients. The number of these
##' additional updates is controlled by \code{refinements} (default one).
##' The final estimates are returned as \code{sigma.x2} and \code{sigma.z2}.
##'
##' As an optional restriction, \code{estimate.inflation = FALSE} fixes both coefficients
##' at one, imposing equal variability across observations and simulations.
##' This bypasses inflation estimation and refinement, so individual forced runs
##' are unnecessary. Polynomial selection and scaling-factor inference remain
##' part of the procedure.
##'
##' The estimated inflation coefficients enter both the final scaling-factor estimate
##' and its estimated covariance \eqn{\widehat\Xi(\widehat f)/N}. Confidence
##' intervals use the resulting asymptotic normal approximation. The residual
##' consistency statistic and its approximate upper-tail normal p-value are
##' returned as \code{diagnostic.stat} and \code{p.value}. Their reference
##' calibration can be affected by estimation of the inflation coefficients and
##' selection of the shrinkage rule. Missing observations are not supported.
##' }
##' @references Li, H. and Li, Y. (2026). Adaptable Fingerprinting with Nonlinear Shrinkage for Climate Change Detection and Attribution under Variance Heterogeneity. Manuscript.
##' @param Xtilde Numeric N by p matrix of ensemble means.
##' @param Y Numeric observation vector of length N.
##' @param mruns Positive integer ensemble sizes, one per forcing.
##' @param ctlruns Numeric m by N matrix of zero-mean control runs.
##'   m may be smaller than N; the sample covariance may be singular.
##' @param ensembles List of p matrices, each n_i by N, containing individual
##'   forced runs. Required only when estimate.inflation is TRUE. Their means
##'   must agree with Xtilde and their row counts with mruns.
##' @param degree Positive integer polynomial degree; default 3, or 1 for linear.
##' @param estimate.inflation Estimate both variance inflation coefficients from data
##'   (default TRUE). FALSE fixes both at one and skips initialization, the
##'   inflation-coefficient search and refinement.
##' @param conf.level Confidence level for marginal normal intervals.
##' @param coefficient.grid Optional matrix with degree columns, each row
##'   containing ell_1,...,ell_k for f(x)=1+sum(ell_j*x^j).
##' @param grid.size Odd integer number of points per coefficient (at least 5), default 21.
##'   The default grid spans the manuscript bounds with additional points near
##'   zero to resolve the bounded admissible region. At most 100000 combinations.
##' @param psi1 Upper spectral domain bound; defaults to 1.05*max(eigen(S)).
##' @param psi2 Upper bound on f(psi1); defaults to max(1.1,0.7*psi1).
##' @param psi3 Minimum polynomial root separation, default 0.001.
##' @param inflation.grid Optional strictly increasing nonnegative grid for sigma_x^2.
##'   Otherwise an adaptive upper bound and a mixed linear/logarithmic grid is used.
##' @param refinements Number of additional inflation/shrinkage updates, default 1.
##' @return
##' An object of classes \code{fingerprintPoly} and \code{fingerprint}, with a \code{plot()} method for
##' standardized residuals.
##' A list containing the following components, where \eqn{N} is the observation
##' dimension and \eqn{p} is the number of forcings:
##' \item{beta}{A named numeric vector of length \eqn{p} of estimated forcing
##' scaling factors \eqn{\widehat\beta}.}
##' \item{var}{The \eqn{p\times p} estimated covariance matrix of
##' \eqn{\widehat\beta}, equal to \eqn{\widehat\Xi(\widehat f)/N}.}
##' \item{sd}{A named numeric vector of standard errors, equal to
##' \code{sqrt(diag(var))}.}
##' \item{ci}{A \eqn{p\times 2} matrix of marginal confidence intervals at
##' \code{conf.level}, with columns \code{lower} and \code{upper}.}
##' \item{sigma.x2}{The variance inflation coefficient for individual forced simulations
##' relative to observations, estimated from data or fixed at one.}
##' \item{sigma.z2}{The variance inflation coefficient for control runs relative to
##' observations, estimated from data or fixed at one.}
##' \item{coefficients}{A numeric vector \eqn{(1,\ell_1,\ldots,\ell_k)} defining
##' the selected polynomial \eqn{\widehat f(t)=1+\sum_{j=1}^k\ell_jt^j}.
##' Coefficients are ordered by ascending power; the intercept is fixed at one.}
##' \item{weight}{The \eqn{N\times N} regularized covariance weight matrix
##' \eqn{W=\widehat f(S)}. Its inverse is used in the weighted TLS loss.}
##' \item{residuals}{A numeric vector of length \eqn{N} of fitted residuals
##' \eqn{Y-\widetilde X\widehat\beta}.}
##' \item{standardized.residuals}{Residuals divided by their estimated marginal
##' standard deviations, computed from
##' \eqn{(1+\widehat\sigma_x^2\widehat\beta^\top D\widehat\beta)
##' S_{ii}/\widehat\sigma_z^2}.}
##' \item{diagnostic.stat}{The standardized residual consistency statistic \eqn{\widehat Q/\sqrt{\widehat V}}, as defined in Li and Li (2026).}
##' \item{p.value}{The approximate upper-tail standard normal p-value for residual
##' consistency. Factor estimation and shrinkage selection can affect its calibration.}
##' \item{call}{The matched function call.}
##' @seealso \code{\link{fingerprint}}, \code{\link{fingerprintLatent}}; \code{\link{plot.fingerprint}} for residual plots;
##' \code{\link{fingerprintLegacy}} for earlier fingerprinting methods.
##' @examples
##' set.seed(12)
##' N <- 60
##' signal <- 3 * sin(seq_len(N) / 6)
##' ensembles <- list(sweep(matrix(rnorm(30 * N), 30), 2, signal, "+"))
##' Xtilde <- matrix(colMeans(ensembles[[1]]), ncol = 1)
##' Z <- matrix(rnorm(120 * N), nrow = 120)
##' Y <- signal + rnorm(N)
##'
##' fit <- fingerprintPoly(Xtilde, Y, mruns = 30, ctlruns = Z,
##'                    ensembles = ensembles)
##'
##' # Observation index, with the default 95% normal reference band.
##' plot(fit)
##'
##' # Use time coordinates and customize the appearance.
##' years <- 1961 + seq_len(N) - 1
##' plot(fit, index = years, xlab = "Year", col = "steelblue",
##'      band.col = "grey90")
##'
##' # Change the pointwise reference level to 90%.
##' plot(fit, reference.level = 0.90)
##' @export
fingerprintPoly <- function(
    Xtilde,
    Y,
    mruns,
    ctlruns,
    ensembles = NULL,
    degree = 3L,
    estimate.inflation = TRUE,
    conf.level = 0.9,
    coefficient.grid = NULL,
    grid.size = 21L,
    psi1 = NULL,
    psi2 = NULL,
    psi3 = 1e-3,
    inflation.grid = NULL,
    refinements = 1L
) {
    Call <- match.call()

    if (!.fp_finite_matrix(Xtilde) || !.fp_finite_matrix(ctlruns))
        stop("Xtilde and ctlruns must be finite numeric matrices")

    N <- nrow(Xtilde)
    p <- ncol(Xtilde)
    m <- nrow(ctlruns)
    if (
        p < 1L || N <= p || m < 2L || ncol(ctlruns) != N ||
        !is.numeric(Y) || length(Y) != N || any(!is.finite(Y))
    )
        stop("Require N > p >= 1, m >= 2 and compatible finite observations/control runs")

    Y <- as.numeric(Y)
    if (
        !is.numeric(mruns) || length(mruns) != p ||
        any(!is.finite(mruns) | mruns < 1 | mruns != floor(mruns))
    )
        stop("mruns must contain p positive integers")
    if (
        !.fp_integer_scalar(degree, 1) ||
        !.fp_integer_scalar(refinements, 0) ||
        !.fp_integer_scalar(grid.size, 5) || grid.size %% 2 != 1
    )
        stop("Invalid degree, refinements or odd grid.size (>= 5)")
    if (
        !is.logical(estimate.inflation) || length(estimate.inflation) != 1 ||
        is.na(estimate.inflation)
    )
        stop("estimate.inflation must be TRUE or FALSE")
    if (!.fp_positive(conf.level) || conf.level >= 1 || !.fp_positive(psi3))
        stop("Require 0 < conf.level < 1 and psi3 > 0")
    if (
        !is.null(inflation.grid) && (
            !is.numeric(inflation.grid) || length(inflation.grid) < 2 ||
            any(!is.finite(inflation.grid)) || any(inflation.grid < 0) ||
            any(diff(inflation.grid) <= 0)
        )
    )
        stop("inflation.grid must be finite, nonnegative and strictly increasing")

    # Control covariance and spectral coordinates.
    D <- diag(1 / mruns, p)
    S <- crossprod(ctlruns) / m
    es <- eigen(S, symmetric = TRUE)
    # Zero eigenvalues are allowed; truncate negative roundoff only.
    ev <- pmax(es$values, 0)
    if (max(ev) <= 0)
        stop("Control covariance must have positive total variation")
    if (is.null(psi1))
        psi1 <- 1.05 * max(ev)
    if (is.null(psi2))
        psi2 <- max(1.1, 0.7 * psi1)
    if (
        !.fp_positive(psi1) || psi1 < max(ev) ||
        !.fp_positive(psi2) || psi2 <= 1
    )
        stop("Require psi1 >= largest sample eigenvalue and psi2 > 1")

    Xt <- crossprod(es$vectors, Xtilde)
    yt <- as.vector(crossprod(es$vectors, Y))

    # Within-ensemble and control variation.
    VX <- NULL
    VZ <- mean(ev)
    if (estimate.inflation) {
        if (!is.list(ensembles) || length(ensembles) != p || sum(mruns) <= p)
            stop("Factor estimation requires individual ensembles and positive within-ensemble degrees of freedom")
        ss <- 0
        for (i in seq_len(p)) {
            E <- ensembles[[i]]
            if (
                !.fp_finite_matrix(E) ||
                !identical(dim(E), c(as.integer(mruns[i]), N))
            )
                stop("Each ensemble must be an mruns[i] by N finite numeric matrix")
            if (!isTRUE(all.equal(
                as.numeric(colMeans(E)),
                as.numeric(Xtilde[, i]),
                tolerance = 1e-07
            )))
                stop("Ensemble means disagree with Xtilde")
            ss <- ss + sum(sweep(E, 2, Xtilde[, i])^2)
        }
        VX <- ss / (N * (sum(mruns) - p))
        if (VX <= 0)
            stop("Within-ensemble variation must be positive")
    }

    # Construct and filter candidate polynomial coefficients.
    if (is.null(coefficient.grid)) {
        if (degree * log(grid.size) > log(1e+05))
            stop("Default grid exceeds 100000 combinations; supply coefficient.grid or reduce grid.size")
        ratio <- (max(ev) + 1) / (min(ev) + 1)
        axes <- lapply(seq_len(degree), function(j) {
            bound <- ratio^(1 - j)
            small <- min(bound / 4, (psi2 - 1) / (2 * degree * psi1^j))
            positive.grid <- exp(seq(
                log(small),
                log(bound),
                length.out = (grid.size - 1) / 2
            ))
            c(-rev(positive.grid), 0, positive.grid)
        })
        coefficient.grid <- as.matrix(expand.grid(axes))
    }
    if (
        !.fp_finite_matrix(coefficient.grid) ||
        ncol(coefficient.grid) != degree || nrow(coefficient.grid) < 1
    )
        stop("coefficient.grid must have degree columns and at least one finite row")
    candidates <- lapply(
        seq_len(nrow(coefficient.grid)),
        function(i) c(1, coefficient.grid[i, ])
    )
    candidates <- candidates[vapply(
        candidates,
        .fp_admissible,
        logical(1),
        degree = degree,
        psi1 = psi1,
        psi2 = psi2,
        psi3 = psi3
    )]
    if (!length(candidates))
        stop("No admissible polynomials; refine coefficient.grid or bounds")

    # Initialize and refine the inflation coefficients and polynomial.
    inflation <- c(x = 1, z = 1)
    if (estimate.inflation) {
        u <- mean(ev)
        denom <- sum((S - diag(u, N))^2)
        numerator <- sum(rowSums(ctlruns^2)^2) - m * sum(S^2)
        q <- if (denom <= 0)
            1
        else min(1, max(0, numerator / (m^2 * denom)))
        if (q <= 0)
            stop("Normalized Ledoit--Wolf initialization has zero intercept")
        a0 <- c(1, (1 - q) / (q * u))
        inflation <- .fp_estimate(
            a0,
            ev = ev,
            Xt = Xt,
            yt = yt,
            N = N,
            mruns = mruns,
            D = D,
            inflation.grid = inflation.grid,
            Y = Y,
            Xtilde = Xtilde,
            VX = VX,
            VZ = VZ
        )
    }
    best <- .fp_select(
        inflation,
        candidates = candidates,
        ev = ev,
        Xt = Xt,
        yt = yt,
        N = N,
        mruns = mruns,
        D = D,
        degree = degree,
        m = m,
        p = p
    )
    if (estimate.inflation && refinements > 0)
        for (i in seq_len(refinements)) {
            inflation <- .fp_estimate(
                best$a,
                ev = ev,
                Xt = Xt,
                yt = yt,
                N = N,
                mruns = mruns,
                D = D,
                inflation.grid = inflation.grid,
                Y = Y,
                Xtilde = Xtilde,
                VX = VX,
                VZ = VZ
            )
            best <- .fp_select(
                inflation,
                candidates = candidates,
                ev = ev,
                Xt = Xt,
                yt = yt,
                N = N,
                mruns = mruns,
                D = D,
                degree = degree,
                m = m,
                p = p
            )
        }

    # Scaling inflation, confidence intervals and residual diagnostics.
    beta <- best$fit$beta
    variance <- best$mom$Xi / N
    se <- sqrt(diag(variance))
    labels <- colnames(Xtilde)
    if (is.null(labels))
        labels <- paste0("forcing", seq_len(p))

    names(beta) <- names(se) <- labels
    dimnames(variance) <- list(labels, labels)
    ci <- cbind(
        lower = beta - stats::qnorm((1 + conf.level) / 2) * se,
        upper = beta + stats::qnorm((1 + conf.level) / 2) * se
    )

    residual <- as.vector(Y - Xtilde %*% beta)
    Q <- sum(crossprod(es$vectors, residual)^2 / best$fit$w) / N - best$mom$scale *
        best$mom$Pi
    V <- 2 / N * best$mom$scale^2 * best$mom$K

    fit <- list(
        beta = beta,
        var = variance,
        sd = se,
        ci = ci,
        sigma.x2 = unname(inflation[1]),
        sigma.z2 = unname(inflation[2]),
        coefficients = best$a,
        weight = tcrossprod(
            sweep(es$vectors, 2, sqrt(best$fit$w), "*")
        ),
        residuals = residual,
        standardized.residuals = residual /
            sqrt(best$mom$scale * diag(S) / inflation[2]),
        diagnostic.stat = Q / sqrt(V),
        p.value = stats::pnorm(Q / sqrt(V), lower.tail = FALSE),
        call = Call
    )
    class(fit) <- c("fingerprintPoly", "fingerprint")
    fit
}

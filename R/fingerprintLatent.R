##' Detection and attribution under a latent factor model
##'
##' Implements the latent-factor optimal fingerprinting framework proposed by
##' Li and Li (2026). The joint errors-in-variables model links the observed climate
##' response, simulated forcing responses, and preindustrial control runs through a
##' common spiked covariance structure. Variance inflation coefficients account for
##' differences in the magnitude of internal variability between observations and
##' simulations. The procedure corrects empirical spike strengths and noise variance,
##' estimates the inflation coefficients, and selects precision weights by minimizing
##' the estimated total asymptotic variance of the forcing scaling-factor estimator.
##'
##' @details
##' \subsection{Joint model for observations and simulations}{
##' Let \eqn{Y\in\mathbb{R}^N} denote the observed climate response and let
##' \eqn{X=(X_1,\ldots,X_p)\in\mathbb{R}^{N\times p}} contain the unknown,
##' deterministic mean responses to \eqn{p} external forcings. For forcing \eqn{i},
##' \eqn{\widetilde X_{ij}} is the response in ensemble member \eqn{j}, and
##' \eqn{n_i} is the ensemble size. The vectors \eqn{Z_1,\ldots,Z_m} are
##' preindustrial control runs. The joint model is
##' \deqn{Y=X\beta+\epsilon,\qquad \epsilon\sim\mathcal N(0,\Sigma),}
##' \deqn{\widetilde X_{ij}=X_i+\eta_{ij},\qquad
##' \eta_{ij}\sim\mathcal N(0,s_x^2\Sigma),\quad
##' j=1,\ldots,n_i,\quad i=1,\ldots,p,}
##' \deqn{Z_j\sim\mathcal N(0,s_z^2\Sigma),\qquad j=1,\ldots,m.}
##' The observation error \eqn{\epsilon}, all forcing-simulation errors
##' \eqn{\eta_{ij}}, and all control runs \eqn{Z_j} are mutually independent.
##' The regression coefficients \eqn{\beta} are the forcing scaling factors used
##' for detection and attribution. The positive variance inflation coefficients
##' \eqn{s_x^2} and \eqn{s_z^2} are variance ratios for forced simulations and
##' controls relative to observations; they need not exceed one. Setting both to
##' one imposes equal internal variability across the three data sources.
##'
##' The observation error has the latent-factor representation
##' \deqn{\epsilon=\sum_{k=1}^{K}\sqrt{\ell_k}U_k f_k+\sigma e,
##' \qquad f_k\stackrel{\mathrm{iid}}{\sim}\mathcal N(0,1),
##' \quad e\sim\mathcal N(0,I_N).}
##' Here \eqn{f_1,\ldots,f_K} are mutually independent latent factors, independent
##' of the idiosyncratic noise vector \eqn{e}. The deterministic vectors
##' \eqn{U_1,\ldots,U_K} are orthonormal loading directions;
##' \eqn{\sqrt{\ell_k}U_k} is the loading vector for factor \eqn{k},
##' \eqn{\ell_k>0} is its strength, and \eqn{\sigma^2>0} is the isotropic
##' idiosyncratic noise variance. Consequently,
##' \deqn{\Sigma=U\operatorname{diag}(\ell_1,\ldots,\ell_K)U^\top+\sigma^2 I_N,
##' \qquad U=(U_1,\ldots,U_K).}
##' The covariance eigenvalue in direction \eqn{U_k} is \eqn{\ell_k+\sigma^2}.
##' Forcing-simulation errors and controls share these loading directions and
##' relative factor strengths, with independent realizations of the latent factors
##' and idiosyncratic noise, scaled by \eqn{s_x} and \eqn{s_z}, respectively.
##'
##' Each unknown forcing response is estimated by its ensemble mean,
##' \deqn{\overline X_i=\frac{1}{n_i}\sum_{j=1}^{n_i}\widetilde X_{ij},
##' \qquad \operatorname{E}(\overline X_i)=X_i,\qquad
##' \operatorname{Cov}(\overline X_i)=\frac{s_x^2}{n_i}\Sigma.}
##' The columns of \code{Xtilde} contain these means, \code{mruns} contains the
##' ensemble sizes, and the rows of \code{ctlruns} contain the controls.
##' Individual forcing ensembles are supplied in \code{ensembles} when estimating
##' the variance inflation coefficients. Uncertainty in the ensemble means gives
##' rise to the errors-in-variables formulation.
##' }
##' \subsection{Scaling-factor estimation}{
##' Write \eqn{\overline X=(\overline X_1,\ldots,\overline X_p)} and
##' \eqn{D=\operatorname{diag}(1/n_1,\ldots,1/n_p)}. For a positive semidefinite
##' precision weight
##' \eqn{W} and a given standard deviation multiplier \eqn{s=s_x}, the scaling
##' factors minimize
##' \deqn{\widehat\beta(s,W)=\arg\min_\beta
##' \frac{(Y-\overline X\beta)^\top W(Y-\overline X\beta)}
##' {1+s^2\beta^\top D\beta}.}
##' The denominator accounts for uncertainty in the simulated ensemble means.
##' The precision weight \eqn{W} enters the numerator directly, unlike the
##' covariance weight returned by \code{fingerprintPoly}. The estimated
##' variance inflation coefficients \eqn{s_x^2} and \eqn{s_z^2} are returned as
##' \code{sigma.x2} and \code{sigma.z2}.
##' }
##' \subsection{Spike correction and weight selection}{
##' Control covariance is \eqn{S=Z^\top Z/m}, without centering. Algorithm 1
##' corrects the control-scale spike strengths \eqn{\pi_k=s_z^2\ell_k} and noise
##' variance \eqn{\sigma_z^2=s_z^2\sigma^2}. The second correction uses the updated
##' noise variance in both the empirical excess eigenvalues and the spike equation.
##' Empirical eigenvectors \eqn{\widetilde U} are retained and their squared
##' alignments with population eigenvectors are estimated from the corrected spikes.
##' The precision weights take the form
##' \deqn{W_\Theta=I_N+\widetilde U\operatorname{diag}(\theta_k)\widetilde U^\top,
##' \qquad \theta_k\geq-1.}
##' A weight of -1 removes that spike direction and is allowed. Equations
##' (10)--(15), including \eqn{-\alpha\Psi_\Theta D} in equation (14), provide the
##' trace and design corrections used in Theorem 4.
##'
##' When \code{n.factors = NULL}, the factor count is estimated from all \eqn{N}
##' ordered eigenvalues of the control covariance, including zero eigenvalues when
##' \eqn{m<N}. For these selectors, the sample size is \eqn{n=m} control runs and
##' the data dimension is \eqn{p=N}; it is not the number of forcing predictors.
##' An explicit nonnegative \code{n.factors} bypasses automatic selection.
##'
##' \code{factor.method = "PLY"} applies the Passemier et al. (2017)
##' criterion. Let \eqn{v_k} be the Algorithm 1 bias-adjusted noise variance when
##' assuming \eqn{k} spikes, and let \eqn{M} be \code{max.factors}. It selects the
##' smallest minimizer of
##' \deqn{v_k+k v_M\frac{p/n+2\sqrt{p/n}}{p}
##' \left(1+\frac{n}{p^{1+\delta}}\right),\quad k=0,\ldots,M,}
##' where \eqn{\delta} is \code{factor.delta}. Noise estimation for trial counts
##' is separated from the eigenvector-alignment checks needed for inference.
##'
##' \code{factor.method = "KN"} (default) applies the sequential Kritchman--Nadler (2008)
##' Tracy--Widom threshold, using \code{RMTstat::qtw(factor.level, beta = 1)}.
##' It tests successive eigenvalues with \eqn{d=p-k}, the finite-sample centering
##' and scale in the supplied selector, and the adjusted noise estimate \eqn{v_k}.
##' It returns the count before the first nonsignificant eigenvalue, or \eqn{M}
##' if all candidates pass. Both methods allow zero factors. The default cap is
##' \eqn{M=\min(m-1,\lfloor 0.1N\rfloor)}; if \eqn{M=0}, selection returns zero
##' without requiring \code{RMTstat}. The KN option otherwise requires that optional
##' package; explicit counts and PLY do not.
##'
##' Selected factors at the detectability boundary have zero alignment, so the
##' inferential variance formula is undefined there. An explicitly smaller
##' \code{n.factors} is needed in that case; the selected count is not silently
##' reduced. Zero factors gives an isotropic model.
##'
##' When \code{theta = NULL}, local L-BFGS-B optimization minimizes
##' \eqn{\operatorname{tr}\widehat\Xi_\Theta}, starting from the naive precision
##' weights and identity weights. This finite numerical procedure does not guarantee
##' a global optimum. Supply \code{theta} for a fixed-weight analysis.
##' }
##' \subsection{Inflation coefficients and inference}{
##' With \code{estimate.inflation = TRUE}, individual forced runs estimate
##' \eqn{\alpha=s_x^2/s_z^2} from the ratio of within-ensemble variation to mean
##' control variation. For each candidate weight, equations (16)--(17) estimate
##' \eqn{s_x} by residual matching and then set \eqn{\widehat s_z=\widehat s_x/\sqrt\alpha}.
##' Both coefficients are re-estimated as the weights change. The API consistently
##' accepts and returns variance multipliers \eqn{s_x^2} and \eqn{s_z^2}.
##' With \code{estimate.inflation = FALSE}, both equal one and individual forced
##' runs are unnecessary. Marginal confidence intervals condition on the selected factor count and use Theorem 4 with covariance
##' \eqn{\widehat\Xi_\Theta/N}.
##'
##' Residual standardization follows equation (21):
##' \eqn{\widehat s_z\widehat\epsilon_i/\sqrt{S_{ii}}}. It omits the extra
##' ensemble-mean variance multiplier used by \code{fingerprintPoly}; the resulting
##' residuals may therefore have extra variability. Plot reference bands are pointwise.
##'
##' With estimated inflation coefficients and exact residual matching, equations (16)--(17) force
##' the residual-test numerator in Section 3.5 to zero. Consequently,
##' \code{diagnostic.stat} is near zero and \code{p.value} near 0.5; they do not
##' provide an independent goodness-of-fit check in this setting. An explicit finite
##' \code{inflation.grid} only approximates this equality. With inflation coefficients fixed at one,
##' the algebraic degeneracy does not apply. Missing values are unsupported.
##' }
##' @references Li, H. and Li, Y. (2026). Fingerprint Analysis for Climate Change Detection and Attribution under a Latent Factor Model. arXiv preprint. \url{https://arxiv.org/abs/2609.13488}.
##' @param Xtilde Numeric \eqn{N\times p} matrix of forcing ensemble means.
##' @param Y Numeric observation vector of length \eqn{N}.
##' @param mruns Positive integer ensemble sizes \eqn{n_i}, one per forcing.
##' @param ctlruns Numeric \eqn{m\times N} matrix of known zero-mean controls. Singular sample covariance is allowed, including \eqn{m<N}.
##' @param n.factors Number \eqn{K} of latent factors. With \code{NULL} (default), estimate it using \code{factor.method}. An explicit nonnegative integer smaller than \eqn{\min(N,m)} overrides selection; zero gives an isotropic model.
##' @param ensembles List of \eqn{p} individual forced ensembles, each \code{mruns[i]} by \eqn{N}. Required to estimate inflation coefficients; their means must agree with \code{Xtilde}.
##' @param estimate.inflation Estimate both variance inflation coefficients (default \code{TRUE}), or fix both at one (\code{FALSE}).
##' @param conf.level Marginal confidence level, default 0.9.
##' @param theta Optional fixed vector of \eqn{K} spike weights, each at least -1. With \code{NULL}, select weights by constrained numerical optimization.
##' @param inflation.grid Optional increasing, nonnegative candidate grid for the variance factor \eqn{s_x^2}, not the standard deviation multiplier \eqn{s_x}. With \code{NULL}, solve residual matching numerically after adaptive bracketing.
##' @param optim.control List of control settings passed to \code{stats::optim} for weight selection.
##' @param factor.method Factor-count selector: \code{"KN"} (default) or \code{"PLY"}. Used when \code{n.factors = NULL}.
##' @param factor.level Tracy--Widom probability for KN selection, default 0.95. Separate from the scaling-factor confidence level \code{conf.level}.
##' @param factor.delta Positive penalty exponent \eqn{\delta} for PLY selection, default 0.05.
##' @param max.factors Largest candidate count. Defaults to \eqn{\min(m-1,\lfloor 0.1N\rfloor)} and must be a nonnegative integer smaller than \eqn{\min(m,N)}. Used only for automatic selection.
##' @return
##' An object of classes \code{fingerprintLatent} and \code{fingerprint},
##' with a \code{plot()} method for standardized residuals. Components are:
##' \item{beta}{Named vector of \eqn{p} forcing scaling-inflation coefficient estimates.}
##' \item{var}{Estimated \eqn{p\times p} coefficient covariance \eqn{\widehat\Xi_\Theta/N}.}
##' \item{sd}{Named vector of standard errors, \code{sqrt(diag(var))}.}
##' \item{ci}{Marginal confidence intervals with columns \code{lower} and \code{upper}.}
##' \item{sigma.x2, sigma.z2}{Estimated variance inflation coefficients \eqn{s_x^2} and
##' \eqn{s_z^2}, or one when \code{estimate.inflation = FALSE}.}
##' \item{n.factors}{The selected or supplied integer factor count.}
##' \item{theta}{Selected or supplied spike weights \eqn{\theta_k}.}
##' \item{spikes}{Bias-corrected control-scale strengths \eqn{\pi_k=s_z^2\ell_k}.}
##' \item{noise.var}{Bias-corrected control noise variance \eqn{s_z^2\sigma^2}.}
##' \item{weight}{\eqn{N\times N} precision weight \eqn{W_\Theta}, used directly
##' in the quadratic loss.}
##' \item{residuals}{Observed minus fitted response, \eqn{Y-\overline X\widehat\beta}.}
##' \item{standardized.residuals}{Residuals standardized according to equation (21).}
##' \item{diagnostic.stat}{Standardized residual consistency statistic from Section 3.5.}
##' \item{p.value}{Upper-tail normal reference p-value; see the estimation limitation in Details.}
##' \item{call}{The matched function call.}
##' @seealso \code{\link{fingerprint}}, \code{\link{fingerprintPoly}}, \code{\link{plot.fingerprint}}
##' @examples
##' set.seed(2026)
##' N <- 80
##' n <- 30
##' m <- 160
##' u <- rep(1 / sqrt(N), N)
##' signal <- 3 * sin(seq_len(N) / 6)
##' Z <- matrix(rnorm(m * N), m) + outer(rnorm(m), sqrt(20) * u)
##' E <- matrix(rnorm(n * N), n) + outer(rnorm(n), sqrt(20) * u)
##' E <- sweep(E, 2, signal, "+")
##' Xtilde <- matrix(colMeans(E), ncol = 1)
##' Y <- signal + rnorm(N) + sqrt(20) * rnorm(1) * u
##' fit <- fingerprintLatent(Xtilde, Y, n, Z,
##'                          ensembles = list(E), factor.method = "PLY")
##' fit$beta
##' fit$ci
##' plot(fit)
##' @export
fingerprintLatent <- function(
    Xtilde,
    Y,
    mruns,
    ctlruns,
    n.factors = NULL,
    ensembles = NULL,
    estimate.inflation = TRUE,
    conf.level = 0.9,
    theta = NULL,
    inflation.grid = NULL,
    optim.control = list(maxit = 200),
    factor.method = c("KN", "PLY"),
    factor.level = 0.95,
    factor.delta = 0.05,
    max.factors = min(nrow(ctlruns) - 1L, floor(ncol(ctlruns) * 0.1))
) {
    Call <- match.call()

    if (!.fp_finite_matrix(Xtilde) || !.fp_finite_matrix(ctlruns)) {
        stop("Xtilde and ctlruns must be finite numeric matrices")
    }
    N <- nrow(Xtilde)
    p <- ncol(Xtilde)
    m <- nrow(ctlruns)
    if (
        p < 1L || N <= p || m < 2L || ncol(ctlruns) != N ||
        !is.numeric(Y) || length(Y) != N || any(!is.finite(Y))
    ) {
        stop("Require N > p >= 1, m >= 2 and compatible finite data")
    }
    Y <- as.numeric(Y)
    if (
        !is.numeric(mruns) || length(mruns) != p ||
        any(!is.finite(mruns) | mruns < 1 | mruns != floor(mruns))
    ) {
        stop("mruns must contain p positive integers")
    }
    factor.method <- match.arg(factor.method)
    if (!is.null(n.factors) && (
        !.fp_integer_scalar(n.factors, 0) || n.factors >= min(N, m)
    )) {
        stop("n.factors must be NULL or a nonnegative integer smaller than min(N, m)")
    }
    if (
        !is.logical(estimate.inflation) || length(estimate.inflation) != 1L ||
        is.na(estimate.inflation)
    ) {
        stop("estimate.inflation must be TRUE or FALSE")
    }
    if (!.fp_positive(conf.level) || conf.level >= 1) {
        stop("conf.level must lie strictly between zero and one")
    }
    if (!is.null(inflation.grid) && (
        !is.numeric(inflation.grid) || length(inflation.grid) < 2L ||
        any(!is.finite(inflation.grid) | inflation.grid < 0) ||
        any(diff(inflation.grid) <= 0)
    )) {
        stop("inflation.grid must be finite, nonnegative and strictly increasing")
    }
    if (!is.list(optim.control)) {
        stop("optim.control must be a list")
    }

    # Control covariance may be singular. Only the selected spike subspace and
    # average trailing eigenvalues are used in the latent model.
    S <- crossprod(ctlruns) / m
    spectrum <- eigen(S, symmetric = TRUE)
    values <- pmax(spectrum$values, 0)
    if (is.null(n.factors)) {
        n.factors <- .lf_select_K(
            values = values,
            n = m,
            p = N,
            method = factor.method,
            conf_level = factor.level,
            delta = factor.delta,
            K_max = max.factors
        )
    }
    if (!is.null(theta) && (
        !is.numeric(theta) || length(theta) != n.factors ||
        any(!is.finite(theta) | theta < -1)
    )) {
        stop("theta must contain n.factors finite weights, each at least -1")
    }

    U <- spectrum$vectors[, seq_len(n.factors), drop = FALSE]
    spike <- .lf_spikes(values, n.factors, m)

    # Estimate alpha = s_x^2 / s_z^2 from individual forcing ensembles.
    alpha <- 1
    if (estimate.inflation) {
        alpha <- .lf_alpha(Xtilde, mruns, ensembles, mean(values))
    }
    design <- list(
        Xtilde = Xtilde, Y = Y, mruns = mruns, U = U,
        N = N, alpha = alpha, spikes = spike$spikes,
        noise = spike$noise, alignment = spike$alignment,
        estimate.inflation = estimate.inflation, inflation.grid = inflation.grid
    )
    best <- .lf_select(design, theta, optim.control)
    if (best$inflation$boundary) {
        warning("Inflation estimate is at inflation.grid's upper endpoint; enlarge the grid")
    }

    # Scaling inflation, intervals, and the paper's residual diagnostic.
    beta <- best$fit$beta
    variance <- best$Xi / N
    se <- sqrt(diag(variance))
    labels <- colnames(Xtilde)
    if (is.null(labels)) labels <- paste0("forcing", seq_len(p))
    names(beta) <- names(se) <- labels
    dimnames(variance) <- list(labels, labels)
    critical <- stats::qnorm((1 + conf.level) / 2)
    ci <- cbind(lower = beta - critical * se, upper = beta + critical * se)
    residual <- best$fit$residuals
    Q <- best$fit$quadratic - best$fit$scale *
        best$trace$delta / best$inflation$z
    V <- 2 / N * best$fit$scale^2 * best$trace$phi / best$inflation$z^2
    precision <- diag(N) + tcrossprod(sweep(U, 2L, best$theta, "*"), U)

    fit <- list(
        beta = beta,
        var = variance,
        sd = se,
        ci = ci,
        sigma.x2 = best$inflation$x,
        sigma.z2 = best$inflation$z,
        n.factors = as.integer(n.factors),
        theta = best$theta,
        spikes = spike$spikes,
        noise.var = spike$noise,
        weight = precision,
        residuals = residual,
        standardized.residuals = sqrt(best$inflation$z) * residual / sqrt(diag(S)),
        diagnostic.stat = Q / sqrt(V),
        p.value = stats::pnorm(Q / sqrt(V), lower.tail = FALSE),
        call = Call
    )
    class(fit) <- c("fingerprintLatent", "fingerprint")
    fit
}

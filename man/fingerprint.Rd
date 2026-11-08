% Documentation for R/fingerprint.R
\name{fingerprint}
\alias{fingerprint}
\title{Detection and attribution of climate change}
\usage{
fingerprint(
    Xtilde,
    Y,
    mruns,
    ctlruns,
    method = c(
        "poly", "linear", "latent", "EE", "PBC", "TS", "TS.TempSt"
    ),
    ...
)
}
\arguments{
\item{Xtilde}{Numeric matrix of forcing ensemble means, one forcing per column.}
\item{Y}{Observation vector or matrix as required by the selected method.}
\item{mruns}{Ensemble sizes, one per forcing.}
\item{ctlruns}{Control runs, with one run per row.}
\item{method}{Estimation method: "poly" (default), "linear", "latent", "EE", "PBC",
  "TS", or "TS.TempSt".}
\item{...}{Arguments for the selected implementation. See fingerprintPoly,
  fingerprintLatent or fingerprintLegacy. In particular, supply ensembles
  when estimating variance inflation. For latent selection, use factor.method
  (PLY or KN) or override the count with n.factors.}
}
\description{
Common entry point for the detection and attribution methods in dacc.
Select polynomial shrinkage under variance heterogeneity, latent-factor
fingerprinting, estimating equations, or two-sample and bootstrap TLS procedures.
}
\details{
The \code{method} argument selects both the scaling-factor estimator and its
uncertainty procedure. Method-specific arguments are passed through \code{...}.
\describe{
\item{\code{"poly"}}{Polynomial-shrinkage fingerprinting proposed by Li and Li
(2026a) in \emph{Adaptable Fingerprinting with Nonlinear Shrinkage for Climate
Change Detection and Attribution under Variance Heterogeneity}.

A polynomial transformation of the control covariance supplies the covariance
weight for weighted total least squares (TLS). Polynomial coefficients minimize
the estimated total asymptotic variance of the scaling-factor estimator;
\code{degree} defaults to 3. Variance inflation coefficients are estimated by
default. Marginal confidence intervals use the estimated asymptotic covariance.}
\item{\code{"linear"}}{Linear-shrinkage fingerprinting obtained as the
degree-one specialization of the method proposed by Li and Li (2026a) in
\emph{Adaptable Fingerprinting with Nonlinear Shrinkage for Climate Change
Detection and Attribution under Variance Heterogeneity}.

A degree-one polynomial supplies the covariance weight for weighted TLS.
Its coefficient minimizes the estimated total asymptotic scaling-factor
variance; any supplied \code{degree} is overridden with 1. Variance inflation
coefficients are estimated by default. Marginal confidence intervals use the
same asymptotic covariance construction as \code{"poly"}. This weight-selection
criterion differs from the covariance-loss criterion of Ledoit and Wolf (2004).}
\item{\code{"latent"}}{Latent-factor fingerprinting proposed by Li and Li
(2026b) in \emph{Fingerprint Analysis for Climate Change Detection and
Attribution under a Latent Factor Model}.

A spiked covariance model represents latent variability plus isotropic noise.
Bias-corrected spike strengths, noise variance and estimated eigenvector
alignments determine a precision weight, selected by minimizing estimated total
asymptotic scaling-factor variance. The factor count uses KN by default
(Kritchman and Nadler, 2008), or PLY (Passemier et al., 2017), unless
\code{n.factors} is supplied. Variance inflation coefficients are estimated by
default. Scaling factors use precision-weighted TLS, and marginal confidence
intervals use the latent model's estimated asymptotic covariance.}
\item{\code{"EE"}}{Estimating-equation fingerprinting developed by Ma et al.
(2023) in \emph{Optimal Fingerprinting with Estimating Equations} and refined
by Li et al. (2025) in \emph{Improved Optimal Fingerprinting Based on Estimating
Equations Reaffirms Anthropogenic Effect on Global Warming}.

The estimating equations correct for measurement error in forcing ensemble
means. A pooled spatial covariance supplies the working weight, with a
block-diagonal temporal-independence structure for multiple time steps.
Marginal confidence intervals use an estimating-equation variance estimate
and a normal approximation. Options include variance calibration through
\code{cal.a}, linear contrasts through \code{Proj}, missing-response handling,
and block-bootstrap variance estimation.}
\item{\code{"PBC"}}{Parametric-bootstrap-calibrated fingerprinting proposed
by Li et al. (2021) in \emph{Uncertainty in Optimal Fingerprinting Is
Underestimated}.

A regularized covariance weight is used for TLS, with linear shrinkage for
\code{cov.method = "l2"} or nonlinear minimum-variance shrinkage for
\code{cov.method = "mv"}. Parametric replicates regenerate the response,
forcing ensemble means and controls from the fitted Gaussian model and
re-estimate the weight. Bootstrap and normal-approximation TLS intervals are
returned together with standard-error calibration ratios in \code{pbc.ratio};
the ratios are supplied separately rather than applied to the interval endpoints.}
\item{\code{"TS"}}{Two-sample regularized fingerprinting described by
Ribes et al. (2013) in \emph{Application of Regularised Optimal Fingerprinting
to Attribution. Part I: Method, Properties and Idealised Analysis}, using the
TLS framework of Allen and Stott (2003).

One control sample estimates the covariance weight used for TLS; a second
sample, \code{ctlruns.bhvar}, supplies the variance estimate. Marginal confidence
intervals use a normal approximation. Nonparametric bootstrap intervals are
also returned, based on resampled coordinates of the transformed regression.
Li et al. (2021) assess the uncertainty properties of the two-sample procedure.}
\item{\code{"TS.TempSt"}}{Two-sample fingerprinting with a temporally
independent working covariance, combining the framework described by Ribes
et al. (2013) in \emph{Application of Regularised Optimal Fingerprinting to
Attribution. Part I: Method, Properties and Idealised Analysis} with covariance
pooling used by Ma et al. (2023) and Li et al. (2025).

Control data are pooled across time to estimate a common spatial covariance,
which is repeated in diagonal temporal blocks and used as the TLS covariance
weight. A second control sample supplies the variance estimate. Marginal
normal-approximation intervals and nonparametric bootstrap intervals follow
the same construction as \code{"TS"}.}
}

For \code{"poly"}, \code{"linear"} and \code{"latent"},
\code{estimate.inflation = FALSE} fixes both variance inflation coefficients at
one; otherwise individual forcing members must be supplied in \code{ensembles}.
The returned \code{weight} is a covariance weight for polynomial/linear fits
and a precision weight for latent fits. The quadratic TLS criterion uses the
inverse of the former and the latter directly. These fits support the common
standardized-residual \code{plot()} method. See the implementation help pages
for the residual diagnostics and their calibration limitations.

The methods \code{"EE"}, \code{"PBC"}, \code{"TS"} and \code{"TS.TempSt"}
are implemented by \code{fingerprintLegacy} with
the selected method passed unchanged. \code{ctlruns} supplies
\code{ctlruns.sigma}; \code{ctlruns.bhvar} defaults to the same sample for
compatibility. Supply a separate second control sample for a two-sample
analysis. Arguments \code{S} and \code{T} describe the spatial and temporal
dimensions used by \code{"EE"} and \code{"TS.TempSt"}. Set \code{B} to a
positive number of bootstrap replicates for \code{"PBC"}, \code{"TS"} and
\code{"TS.TempSt"}. The modern options \code{degree}, \code{estimate.inflation}
and latent-factor selection controls do not apply to these four methods.
}
\value{
The result returned by the selected implementation. Polynomial and
  latent fits inherit from class fingerprint and support plot(). The EE and
  TLS procedures retain their implementation-specific list structures.
}
\references{
Li, H. and Li, Y. (2026a). Adaptable Fingerprinting with Nonlinear Shrinkage
for Climate Change Detection and Attribution under Variance Heterogeneity.
arXiv preprint. \url{https://arxiv.org/abs/2608.12496}.

Li, H. and Li, Y. (2026b). Fingerprint Analysis for Climate Change Detection
and Attribution under a Latent Factor Model. arXiv preprint. \url{https://arxiv.org/abs/2609.13488}.

Ma, S., Wang, T., Yan, J. and Zhang, X. (2023). Optimal Fingerprinting with
Estimating Equations. \emph{Journal of Climate}, \bold{36}, 7109--7122.
\doi{10.1175/JCLI-D-22-0681.1}.

Li, Y., Wang, T., Yan, J. and Zhang, X. (2025). Improved Optimal Fingerprinting
Based on Estimating Equations Reaffirms Anthropogenic Effect on Global Warming.
\emph{Journal of Climate}, \bold{38}, 1779--1790.
\doi{10.1175/JCLI-D-24-0193.1}.

Li, Y., Chen, K., Yan, J. and Zhang, X. (2021). Uncertainty in Optimal
Fingerprinting Is Underestimated. \emph{Environmental Research Letters},
\bold{16}, 084043.

Allen, M. R. and Stott, P. A. (2003). Estimating Signal Amplitudes in Optimal
Fingerprinting, Part I: Theory. \emph{Climate Dynamics}, \bold{21}, 477--491.
\doi{10.1007/s00382-003-0313-9}.

Ribes, A., Planton, S. and Terray, L. (2013). Application of Regularised
Optimal Fingerprinting to Attribution. Part I: Method, Properties and Idealised
Analysis. \emph{Climate Dynamics}, \bold{41}, 2817--2836.
\doi{10.1007/s00382-013-1735-7}.

Ledoit, O. and Wolf, M. (2004). A Well-Conditioned Estimator for
Large-Dimensional Covariance Matrices. \emph{Journal of Multivariate Analysis},
\bold{88}, 365--411. \doi{10.1016/S0047-259X(03)00096-4}.

Kritchman, S. and Nadler, B. (2008). Determining the Number of Components in a
Factor Model from Limited Noisy Data. \emph{Chemometrics and Intelligent
Laboratory Systems}, \bold{94}, 19--32.

Passemier, D., Li, Z. and Yao, J. (2017). On Estimation of the Noise Variance
in High Dimensional Probabilistic Principal Component Analysis.
\emph{Journal of the Royal Statistical Society: Series B}, \bold{79}, 51--67.
}
\examples{
set.seed(12)
N <- 60
X <- matrix(3 * sin(seq_len(N) / 6), ncol = 1)
Z <- matrix(rnorm(120 * N), nrow = 120)
Y <- as.vector(X + rnorm(N))
fit <- fingerprint(X, Y, mruns = 30, ctlruns = Z,
                   method = "linear", estimate.inflation = FALSE)
fit$ci
plot(fit)
}
\seealso{\code{\link{fingerprintPoly}}, \code{\link{fingerprintLatent}},
\code{\link{fingerprintLegacy}}, \code{\link{plot.fingerprint}}}

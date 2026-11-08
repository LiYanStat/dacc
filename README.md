# dacc: Detection and Attribution Analysis of Climate Change

`dacc` provides a unified statistical framework for detection and attribution
(D&A) analysis of climate change. The package estimates scaling factors that
relate observed climate responses to simulated responses under external forcings
and quantifies their uncertainty. Supported approaches include regularized total
least squares, polynomial covariance shrinkage, latent-factor modeling,
estimating equations, and bootstrap calibration. Complementary functions support
preprocessing of gridded climate data and regularized covariance estimation.

## Installation

Install the development version from a local source directory:

```r
install.packages("remotes") # if needed
remotes::install_local("path/to/dacc", dependencies = TRUE)
library(dacc)
```

The default latent-factor selection procedure, KN, requires the optional package
`RMTstat`. It can be installed with `install.packages("RMTstat")`. PLY selection
and analyses with a specified factor count do not require this dependency.

## Statistical methods

The primary interface is `fingerprint()`. Its `method` argument selects the
estimation and uncertainty-quantification procedure; the default is `"poly"`.

| Method | Estimation and uncertainty | Main references |
|---|---|---|
| `"poly"` | Weighted TLS with polynomial covariance shrinkage selected by estimated total asymptotic scaling-factor variance; asymptotic confidence intervals. Degree 3 by default. | Li and Li (2026a) |
| `"linear"` | Degree-one specialization of polynomial fingerprinting, with the same variance criterion and asymptotic inference. Always sets degree to 1. | Li and Li (2026a) |
| `"latent"` | Precision-weighted TLS under a latent-factor covariance model, with bias-corrected spikes, adaptive weights, and asymptotic confidence intervals. | Li and Li (2026b) |
| `"EE"` | Measurement-error-corrected estimating equations with a pooled spatial working covariance and normal-approximation confidence intervals. | Ma et al. (2023); Li et al. (2025) |
| `"PBC"` | Regularized TLS with parametric bootstrap calibration; returns TLS intervals and separate standard-error calibration ratios. | Li et al. (2021) |
| `"TS"` | Two-sample regularized TLS, using separate control samples for weighting and uncertainty estimation; normal and nonparametric bootstrap intervals. | Allen and Stott (2003); Ribes et al. (2013) |
| `"TS.TempSt"` | Two-sample TLS with a pooled spatial covariance repeated across temporal blocks under working temporal independence. | Ribes et al. (2013); Ma et al. (2023); Li et al. (2025) |

Method-specific interfaces are available through `fingerprintPoly()` and
`fingerprintLatent()`. The estimating-equation and bootstrap procedures are
implemented in `fingerprintLegacy()` and selected directly through
`fingerprint()`. The R documentation (`?fingerprint`) describes each procedure
and provides links to its implementation and methodological references.

## Data structure

Fingerprinting requires three sources of information: observed climate
responses, ensembles of simulations under specified external forcings, and
preindustrial control simulations representing internal climate variability.

| Input | Contents | Shape |
|---|---|---|
| `Y` | Observed climate response, stacked over the selected spatial and temporal coordinates. | Length `N`; supply an `N x 1` matrix for the TLS bootstrap procedures. |
| `Xtilde` | Ensemble mean response to each of `p` external forcings. | `N x p` matrix. |
| `mruns` | Number of ensemble members contributing to each forcing mean. | Length `p`. |
| `ctlruns` | Internal-variability control runs, one per row. | `m x N` matrix. |
| `ensembles` | Individual forced simulations, needed when estimating variance inflation coefficients. | List of `p` matrices, each `mruns[i] x N`. |

All data must be defined on the same spatial and temporal coordinates and use
consistent ordering. When individual forcing ensembles are supplied, their
means must agree with `Xtilde`. The polynomial and latent-factor procedures
assume zero-mean control simulations and accommodate singular sample covariance
matrices, including settings with fewer control runs than observation
coordinates (`m < N`).

The underlying regression model is

$$
Y=X\beta+\epsilon,
$$

where the columns of $X$ are the unknown mean responses to individual external
forcings and $\beta$ contains their scaling factors. Ensemble means in `Xtilde`
estimate these responses. Their sampling uncertainty gives rise to an
errors-in-variables regression, while the covariance of $\epsilon$ represents
dependence across observation coordinates. Confidence intervals for the scaling
factors provide the basis for detection and attribution statements.

The `"poly"`, `"linear"`, and `"latent"` methods allow simulation variability to
differ from observational variability through the variance inflation coefficients
$s_x^2$ and $s_z^2$, for forced simulations and controls, respectively. Both are
estimated by default. Setting `estimate.inflation = FALSE` imposes the
restriction $s_x^2=s_z^2=1$; individual forcing ensembles are then unnecessary.
These coefficients are distinct from the latent factors used to represent the
covariance structure in the `"latent"` method.

## Documentation

Method-specific assumptions, numerical controls, return values, and diagnostic
limitations are documented in `?fingerprint`, `?fingerprintPoly`,
`?fingerprintLatent`, and `?fingerprintLegacy`. See `?fpPrep` for gridded-data
preparation, `?Covest` for covariance estimation, and `?simDat` and `?globalDat`
for the included datasets.

## Example analysis

The following example simulates a single-forcing analysis with a dominant
latent variability component and fits the polynomial-shrinkage model.

```r
set.seed(2026)
N <- 80
n <- 30
m <- 160
u <- rep(1 / sqrt(N), N)
signal <- 3 * sin(seq_len(N) / 6)

Z <- matrix(rnorm(m * N), m) + outer(rnorm(m), sqrt(20) * u)
E <- matrix(rnorm(n * N), n) + outer(rnorm(n), sqrt(20) * u)
E <- sweep(E, 2, signal, "+")
Xtilde <- matrix(colMeans(E), ncol = 1,
                 dimnames = list(NULL, "forcing"))
Y <- signal + rnorm(N) + sqrt(20) * rnorm(1) * u

fit <- fingerprint(Xtilde, Y, mruns = n, ctlruns = Z,
                   method = "poly", ensembles = list(E))
fit$beta
fit$ci
c(fit$sigma.x2, fit$sigma.z2)
plot(fit)
```

Polynomial and latent-factor fits contain scaling-factor estimates, their
estimated covariance matrix and standard errors, marginal confidence intervals,
variance inflation coefficient estimates, and residual diagnostics. The S3
method `plot()` displays standardized residuals with a reference line at zero
and a pointwise standard normal reference band. Complete return-value
definitions are provided in the method-specific help pages.

Alternative weighting models can be fitted to the same data:

```r
linear <- fingerprint(Xtilde, Y, n, Z, method = "linear", ensembles = list(E))

# PLY estimates the factor count without an optional dependency.
latent <- fingerprint(Xtilde, Y, n, Z, method = "latent",
                      ensembles = list(E), factor.method = "PLY")
latent$n.factors

# Fix both variance inflation coefficients at one.
linear.equal <- fingerprint(Xtilde, Y, n, Z, method = "linear",
                            estimate.inflation = FALSE)
```

## References

- Allen, M. R. and Stott, P. A. (2003). *Estimating Signal Amplitudes in Optimal
  Fingerprinting, Part I: Theory*. Climate Dynamics, 21, 477–491.
- Ribes, A., Planton, S. and Terray, L. (2013). *Application of Regularised
  Optimal Fingerprinting to Attribution. Part I: Method, Properties and Idealised
  Analysis*. Climate Dynamics, 41, 2817–2836.
- Ma, S., Wang, T., Yan, J. and Zhang, X. (2023). *Optimal Fingerprinting with
  Estimating Equations*. Journal of Climate, 36, 7109–7122.
- Li, Y., Wang, T., Yan, J. and Zhang, X. (2025). *Improved Optimal Fingerprinting
  Based on Estimating Equations Reaffirms Anthropogenic Effect on Global Warming*.
  Journal of Climate, 38, 1779–1790.
- Li, Y., Chen, K., Yan, J. and Zhang, X. (2021). *Uncertainty in Optimal
  Fingerprinting Is Underestimated*. Environmental Research Letters, 16, 084043.
- Li, H. and Li, Y. (2026a). *Adaptable Fingerprinting with Nonlinear Shrinkage
  for Climate Change Detection and Attribution under Variance Heterogeneity*.
  [arXiv:2608.12496](https://arxiv.org/abs/2608.12496).
- Li, H. and Li, Y. (2026b). *Fingerprint Analysis for Climate Change Detection
  and Attribution under a Latent Factor Model*. [arXiv:2609.13488](https://arxiv.org/abs/2609.13488).
- Kritchman, S. and Nadler, B. (2008). *Determining the Number of Components in
  a Factor Model from Limited Noisy Data*. Chemometrics and Intelligent
  Laboratory Systems, 94, 19–32.
- Passemier, D., Li, Z. and Yao, J. (2017). *On Estimation of the Noise Variance
  in High Dimensional Probabilistic Principal Component Analysis*. Journal of
  the Royal Statistical Society: Series B, 79, 51–67.

Additional methodological references are listed in the corresponding R help
pages.

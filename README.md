# dacc: Detection and attribution of climate change

Implements the method proposed in work by Li and Li (2026) for detection
and attribution (D&A) analysis of climate change using polynomial shrinkage
under variance heterogeneity. The framework allows forced simulations and
control runs to have different internal variability magnitudes from observations,
while sharing a common covariance structure. It estimates forcing scaling factors,
optionally estimates variance inflation factors, and provides asymptotic confidence
intervals and approximate residual diagnostics.

## Function roles

| Function | Role |
|---|---|
| `fingerprint()` | Implements the method proposed in work by Li and Li (2026): D&A analysis using polynomial shrinkage under variance heterogeneity. Defaults to degree 3 and estimated inflation factors. |
| `fingerprintLegacy()` | Earlier EE, PBC, TS and TS.TempSt methods, retained for reproducing previous analyses and comparing methods. Preserves the original estimation procedures, confidence intervals, arguments and defaults. |

These are separate functions. Setting `degree = 1` and
`estimate.factors = FALSE` in `fingerprint()` keeps the new shrinkage selection
and inference procedure; it does not invoke an earlier method.

## Model and framework

For observations $Y$, true forcing responses $X_i$, individual forced simulations
$\widetilde X_{ik}$, and control runs $Z_j$, Li and Li (2026) consider

$$
Y=X\beta+\epsilon,\quad \epsilon\sim\mathcal N(0,\Sigma),
$$
$$
\widetilde X_{ik}=X_i+\eta_{ik},\quad
\eta_{ik}\sim\mathcal N(0,\sigma_x^2\Sigma),\qquad
Z_j\sim\mathcal N(0,\sigma_z^2\Sigma).
$$

Errors are mutually independent. The inflation factors describe forced-simulation
and control variability relative to observations while preserving a common
covariance structure. Setting both factors to one gives the conventional equal
variability model. Ensemble means form $\widetilde X$, with
$D={diag}(1/n_1,\ldots,1/n_p)$.

The method constructs $W=f(S)$ from $S=m^{-1}\sum_j Z_jZ_j^\top$, using an
admissible polynomial $f(t)=1+\sum_{j=1}^k\ell_jt^j$, and estimates

$$
\widehat\beta(s,W)=\arg\min_\beta
\frac{(Y-\widetilde X\beta)^\top W^{-1}(Y-\widetilde X\beta)}
{1+s\beta^\top D\beta},\qquad s=\sigma_x^2.
$$

Polynomial coefficients minimize estimated total asymptotic variance over a finite
grid. Degree 1 gives linear shrinkage; degree 3 is the recommended setting.
Inflation factors can be estimated from within-ensemble variation and residuals,
or both fixed at one. Confidence intervals use the estimated asymptotic covariance;
the residual consistency diagnostic has an approximate calibration. See
`?fingerprint` for the factor-estimation formulas and numerical controls.

## Manuscript method

```r
# Xtilde: N x p ensemble means; Y: length N; Z: m x N zero-mean controls.
# ensembles: list of p matrices, each n_i x N, with means matching Xtilde.
fit <- fingerprint(Xtilde, Y, mruns, ctlruns = Z,
                       ensembles = ensembles,
                       degree = 3, estimate.factors = TRUE)

# Equal variability: no individual ensemble members or factor search needed.
traditional <- fingerprint(Xtilde, Y, mruns, ctlruns = Z,
                               degree = 1, estimate.factors = FALSE)

# Polynomial degree and factor estimation are independent choices.
cubic.equal <- fingerprint(Xtilde, Y, mruns, ctlruns = Z,
                               degree = 3, estimate.factors = FALSE)
linear.estimated <- fingerprint(Xtilde, Y, mruns, ctlruns = Z,
                                    ensembles = ensembles,
                                    degree = 1, estimate.factors = TRUE)
fit$beta
fit$ci
c(fit$sigma.x2, fit$sigma.z2)
plot(fit)
```

The linear case uses the traditional equal-variability TLS model when factors
are fixed, but selects the linear shrinkage coefficient by the manuscript's
estimated total variance criterion. It need not reproduce `Covest(method="l2")`.
The default grid uses 21 points per coefficient (9,261 combinations for degree 3
before admissibility filtering). The finite grid approximates the manuscript
optimization; increase
`grid.size` or supply `coefficient.grid` to assess search resolution. Higher
polynomial degrees generally require an explicit grid to avoid exponential cost.
Control covariance is `crossprod(Z)/nrow(Z)`, using the model's known zero mean.
Confidence intervals use the manuscript asymptotic covariance, divided by N.
The residual diagnostic p-value is approximate, especially with estimated factors.

## Original fingerprinting methods

Use `fingerprintLegacy()` for the original EE, PBC, TS, and TS.TempSt methods,
with their original arguments and confidence interval procedures:

```r
legacy <- fingerprintLegacy(Xtilde, Y, mruns,
                            ctlruns.sigma = Z, ctlruns.bhvar = Z,
                            S = S, T = T, method = "EE")
```

Existing calls to the original `fingerprint()` should be renamed to
`fingerprintLegacy()`. The new `fingerprint()` uses `ctlruns` for the control
matrix and does not take `S`, `T`, or a legacy `method` argument.

## Reference

Li, H. and Li, Y. (2026). Adaptable Fingerprinting with Nonlinear Shrinkage for Climate Change Detection and Attribution under Variance Heterogeneity. Manuscript.

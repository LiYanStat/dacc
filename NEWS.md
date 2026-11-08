# dacc 0.2.0

## Method interface

* Renamed the polynomial implementation to `fingerprintPoly()`. `fingerprint()`
  is now a common entry point with `method = "poly"` (default), `"linear"`, `"latent"`,
  `"EE"`, `"PBC"`, `"TS"`, or `"TS.TempSt"`. Method-specific arguments are
  passed through `...`.
* Retained `fingerprintLegacy()` and its original procedures. These can be
  selected directly by `method` in `fingerprint()`. Removed the redundant
  `"legacy"` choice and the `legacy.method` selector.

* Added `method = "linear"` as an explicit degree-one polynomial method,
  retaining both estimated and fixed inflation-coefficient options. Any supplied
  degree is overridden with one.

## Latent factor framework

* Added `fingerprintLatent()` implementing the latent-factor framework of
  Li and Li (2026), including Algorithm 1 spike and noise corrections,
  precision-weight selection, variance inflation estimation, asymptotic
  confidence intervals, and residual diagnostics.
* Added automatic factor-count selection from control-covariance eigenvalues
  using Kritchman and Nadler (2008), `factor.method = "KN"` (default), or
  Passemier et al. (2017), `factor.method = "PLY"`. KN uses optional package
  RMTstat.
  An explicit `n.factors` overrides selection. The selected count is returned
  as `n.factors`. Both selectors allow zero factors and configurable search caps.
* Supports singular control covariance and theta weights equal to -1.
* Both modern fit classes inherit from `fingerprint` and use the existing plot
  method, including the pointwise normal residual reference band.
* Documented that exact residual matching with estimated latent inflation
  factors makes the plug-in residual diagnostic zero (p-value near 0.5).
* Added independent matrix checks for trace approximations, weighted TLS,
  covariance estimates, and method dispatch.

## Terminology

* Use variance inflation coefficients for the variability quantities, reserving
  factor terminology for the latent structure. The preferred arguments are now
  `estimate.inflation` and `inflation.grid`; the old argument names are no
  longer accepted. Latent-selection controls (`n.factors`,
  `factor.method`, `factor.level`, `factor.delta`, and `max.factors`) are unchanged.

# dacc 0.1.0

## New features

* Added the method proposed by Li and Li (2026) for detection and attribution
  analysis using polynomial shrinkage under variance heterogeneity, available
  through `fingerprint()`.
* Added estimation of variance inflation factors for forced simulations and
  control runs relative to observations. Set `estimate.factors = FALSE` to fix
  both factors at one and skip factor estimation and refinement.
* Added configurable polynomial degree (`degree`, default 3; degree 1 gives
  linear shrinkage) and coefficient search controls. The default grid uses 21
  points per coefficient, giving 9,261 combinations for degree 3 before filtering.
* Added asymptotic confidence intervals, standardized residuals and an approximate
  residual consistency diagnostic. The diagnostic outputs
  `diagnostic.stat` and `p.value` are returned directly.

* Added a `plot()` method for `fingerprint` fits to display standardized
  residuals against observation index or supplied coordinates.

## Interface changes

* Renamed the original `fingerprint()` to `fingerprintLegacy()`, preserving its
  EE, PBC, TS and TS.TempSt methods, arguments and defaults. Existing analyses
  using these methods should call `fingerprintLegacy()`.
* The new `fingerprint()` uses `ctlruns` for controls and does not require `S`
  or `T`. Individual forced ensemble members are required when estimating
  variance inflation factors.

## Documentation and maintenance

* Documented the model, estimation procedure, arguments and returned components,
  with attribution to Li and Li (2026).
* Added Haoran Li as an author.
* Moved numerical helpers to a separate internal script and added numerical
  regression checks, including cases with singular control covariance.

# dacc 0.0-7

* Added an alternative two sample approach using the weight matrix from EE method

* Fixed error in R code

# dacc 0.0-6

* Added Projection matrix for multiple way analysis with combined forcings

# dacc 0.0-5

* Fixed errors

# dacc 0.0-4

* Added a dataset for global mean temperature

* Updated the main function fingerprint() to include different methods and fix name of inputs

# dacc 0.0-3

## Minor changes

* Added reference for the package

# dacc 0.0-2

## Minor changes

* Fixed reference for functions.


# dacc 0.0-1

## New features

* The first version of package that provides functions for detection and 
  attribution analysis of climate change.

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

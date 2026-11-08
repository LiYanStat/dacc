## Update

This updates dacc from CRAN version 0.0-7 to version 0.2.0. The update adds
polynomial shrinkage and latent factor methods, a unified method interface,
and residual plots. The original procedures remain available through explicit
method choices and fingerprintLegacy(). See NEWS.md for interface changes.

Yan Li remains the designated maintainer.

## Verified check environment

- Local macOS Tahoe 26.6, Apple Silicon, R 4.6.1 (2026-06-24).
- Fresh source build followed by `R CMD check --as-cran dacc_0.2.0.tar.gz`.
- Result: 0 errors, 0 warnings, 1 NOTE.
- Installation, examples, all five test scripts, and PDF manual creation passed.
- CRAN incoming feasibility and DESCRIPTION checks passed.

The NOTE concerns local HTML checking tools: HTML Tidy is too old and the
optional V8 package is unavailable, so HTML validation and mathematical
rendering checks were skipped. These checks remain to be completed.

## Reverse dependencies

The CRAN source package index retrieved on 2026-11-08 lists no reverse
dependencies, including suggested or enhanced dependencies.

## Outstanding before submission

- Obtain current Windows/Linux and R-devel check results. No such results
  have been verified for this release; older claims have been removed.
  A four-job GitHub Actions workflow is prepared in
  .github/workflows/R-CMD-check.yaml but has not yet run.
- Complete HTML validation and mathematical rendering checks.
- Add a public URL for the Li and Li (2026) manuscripts cited in DESCRIPTION,
  as required by CRAN policy for unpublished references.
- Have the designated maintainer review and approve the release and submission.

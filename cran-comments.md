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

## Verified GitHub Actions results

Run: https://github.com/LiYanStat/dacc/actions/runs/42532617190
Checked commit: c8cda2e4e01b1f3cf5345f36567611204fca354f (2026-11-08).

All four jobs completed with Status: OK (0 errors, 0 warnings, 0 NOTEs):

- Ubuntu 26.04.1 LTS, R 4.6.2.
- Ubuntu 26.04.1 LTS, R-devel (2026-11-07 r90850).
- Windows Server 2022 x64, R 4.6.2 UCRT.
- Windows Server 2022 x64, R-devel (2026-11-07 r90850 UCRT).

These jobs used --as-cran --no-manual and disabled CRAN incoming network
checks through the r-lib check action defaults. They do not establish that
HTML validation or PDF manual checks ran. The local full check above covers
PDF manual generation and incoming feasibility. Manuscript URLs were added
after this checked commit and are not covered by these results.

## Reverse dependencies

The CRAN source package index retrieved on 2026-11-08 lists no reverse
dependencies, including suggested or enhanced dependencies.

## Outstanding before submission

- Complete HTML validation and mathematical rendering checks.
- Recheck the final submission tarball after the manuscript URL updates.
- Have the designated maintainer review and approve the release and submission.

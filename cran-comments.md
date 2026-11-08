## Update

This updates dacc from CRAN version 0.0-7 to version 0.2.0. The update adds
polynomial shrinkage and latent factor methods, a unified method interface,
and residual plots. The original procedures remain available through explicit
method choices and fingerprintLegacy(). See NEWS.md for interface changes.

Yan Li remains the designated maintainer.

## Verified check environment

- Local macOS Tahoe 26.6, Apple Silicon, R 4.6.1 (2026-06-24).
- Fresh source build followed by `R CMD check --as-cran dacc_0.2.0.tar.gz`.
- Result: Status: OK (0 errors, 0 warnings, 0 NOTEs).
- Installation, examples, all five test scripts, and PDF/HTML manual checks passed.
- CRAN incoming feasibility and DESCRIPTION checks passed.

HTML validation and mathematical rendering were enabled using temporary
HTML Tidy 5.8.0 and V8 8.2.1 installations. No manual checks were skipped.
The full check was repeated after converting the DESCRIPTION arXiv citations
to their arXiv DOI identifiers; incoming feasibility then passed without a NOTE.

## Verified GitHub Actions results

Run: https://github.com/LiYanStat/dacc/actions/runs/42533622340
Checked commit: 3c2a62d9561e1e1aa5cb5d8cf6efef43c2207838 (2026-11-08).

All four jobs completed with Status: OK (0 errors, 0 warnings, 0 NOTEs):

- Ubuntu 26.04.1 LTS, R 4.6.2.
- Ubuntu 26.04.1 LTS, R-devel (2026-11-07 r90850).
- Windows Server 2022 x64, R 4.6.2 UCRT.
- Windows Server 2022 x64, R-devel (2026-11-07 r90850 UCRT).

These jobs used --as-cran --no-manual and disabled CRAN incoming network
checks through the r-lib check action defaults. They do not establish that
HTML validation or PDF manual checks ran. The local full check above covers
PDF and HTML manual validation and incoming feasibility. The subsequent
DESCRIPTION-only conversion from arXiv URLs to arXiv DOI identifiers is covered
by the local full check, rather than these GitHub results.

## Reverse dependencies

The CRAN source package index retrieved on 2026-11-08 lists no reverse
dependencies, including suggested or enhanced dependencies.

## Outstanding before submission

- Have the designated maintainer review and approve the release and submission.

## CRAN Submission Comments

This is a re-submission of **dacc**.

The previous version was archived on CRAN on 2025-09-08 because reported
errors were not corrected in time. This update fixes those issues and now
passes all CRAN checks with 0 ERRORs and 0 WARNINGs.

### Changes since archival

- Fixed test failures on R-devel that caused the package to be archived.
- Standardized `DESCRIPTION` text and citation formatting.
- Removed the `Date` field to avoid the "Date field is over a month old" NOTE.
- Added `inst/WORDLIST` to whitelist proper names ("Ledoit", "et", "al")
  and avoid spell-check NOTEs.
- Verified examples run within CRAN time limits.

### Check environments

- Local macOS (R 4.4.1): `R CMD check --as-cran` — **OK**
- win-builder (R-devel, R-release): **OK**
- GitHub Actions CI (Ubuntu, Windows, macOS; R-oldrel, R-release, R-devel): **OK**

There are no remaining NOTEs except the expected message that this is a
"new submission" for an archived package.
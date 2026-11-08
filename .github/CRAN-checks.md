# Cross-platform CRAN checks

The R-CMD-check workflow checks Windows and Linux with both R-release and
R-devel. It runs on pushes to main, master, and fingerprint-0.1.0, on pull
requests, and through GitHub's manual workflow dispatch.

To run the prepared checks, push the workflow and current package source to
the repository. Open the repository's Actions tab, select R-CMD-check, and
inspect all four jobs. For subsequent manual runs, select Run workflow and
the branch to check; GitHub exposes manual dispatch after the workflow is
present on the default branch. Push-triggered runs work on the listed branches.

Each job installs all package dependencies, including RMTstat, and runs
R CMD check with --as-cran --no-manual. Errors and warnings fail the job.
NOTEs still require review even when the job succeeds. Check-result artifacts
are retained for successful runs and for failures during the check step.

Linux jobs install HTML Tidy, and all jobs install V8, to enable HTML and
mathematical rendering checks where supported by the R version. Confirm in
the check logs that these checks actually ran rather than being skipped.

PDF manual generation is omitted from these jobs; it passed in the local
full --as-cran check. The r-lib check action disables CRAN incoming network
checks by default. Incoming feasibility passed locally, and should be rerun
on the final submission tarball.

Before recording results in cran-comments.md, download the artifacts and
record the checked commit, actual R and operating-system versions, and the
errors, warnings, and NOTEs from each 00check.log. A green Actions run alone
does not establish that all CRAN requirements are satisfied.

No checks have been run on GitHub yet for this workflow. This setup does not
submit the package to CRAN or upload it to win-builder.

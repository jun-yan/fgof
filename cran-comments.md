# Draft for fgof 0.2-2: not submitted

This is a revival of fgof 0.2-1, archived on 2022-04-27 because check issues
were not corrected in time. The archival example failure is resolved by
normalizing qgrid before scalar conditions; simulation's analogous default
failure and missing namespace imports are also repaired.

Additional regression-tested repairs preserve named/fixed distribution
parameters, evaluate analytical CDF derivatives at the requested points,
restore the missing multivariate empirical CDF without a dependency, and align
bootstrap statistic labels. The public signature and return layout are
unchanged. No runtime dependencies or statistical procedures were added.
The regression suite uses testthat as a suggested test dependency.
The maintainer and contact address are unchanged.

At the maintainer's request, p-values now use (count + 0.5) / (nsim + 1), with
count including ties. This finite-simulation convention avoids exact zero and
one. It changes only the reported p-values relative to the archived formula;
estimates, statistics, and bootstrap draws are unchanged on compatible paths.

## Local checks

- R 4.6.1, x86_64 macOS Sonoma 14.8.7.
- Build succeeded. The testthat suite, namespace-isolation script, and all
  examples passed; current infrastructure check results are saved separately
  from the earlier base R test logs.
- Supplementary R CMD check --as-cran completed with status OK, including PDF
  and HTML manuals, with incoming checks and remote clock verification disabled.
- Full incoming feasibility could not run because R could not access CRAN
  package indexes. This limitation has not been treated as a successful check.
- The pushed pre-adjustment commit 412f5c4 passed all five GitHub Actions jobs:
  Windows/macOS release, Linux release/devel/oldrel-1. The p-value adjustment
  still needs its own pushed CI run.
- A second supplementary check using current CRAN mvtnorm 1.4.2 also passed.
  The current CRAN index listed no direct reverse dependencies in Depends,
  Imports, LinkingTo, Suggests, or Enhances.
- After the approved p-value adjustment, all 12 testthat cases and 111
  expectations passed. The supplementary package check again completed with
  status OK, including tests, examples, and PDF/HTML manuals. Logs are saved
  under misc/audit/pvalue-adjustment/.

## Before using this draft

- Complete normal online --as-cran checks with current R and dependencies, and
  cross-platform checks on the final p-value-adjustment commit.
- Refresh reverse-dependency inspection and these results before submission.
- Obtain the maintainer's explicit approval to submit. No submission has occurred.

The 200-dataset-per-model simulation is a pilot with reported Monte Carlo
intervals, not a claim of exact finite-sample calibration.

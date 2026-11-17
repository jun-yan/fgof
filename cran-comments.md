# Draft for fgof 0.2-2: not submitted

This revives fgof 0.2-1, archived on 2022-04-27 after check issues were not
corrected in time. The qgrid and analogous simulation default failures are
resolved by normalizing choices before scalar conditions. Namespace imports,
fixed/named parameter forwarding, callback evaluation points, multivariate
empirical CDF counts, and bootstrap statistic labels are repaired.

Default multiplier-bootstrap numerical derivatives now use generic scale,
feasibility, stability, and representability checks for standard and custom
models. Explicit numerical controls and analytical callbacks are preserved.
The public signature and returned structure, optimizer, covariance/influence
formulas, and bootstrap procedures are unchanged. No runtime dependencies were
added. Unreliable numerical derivatives can now report an actionable error.

At the maintainer's request, p-values use (count + 0.5) / (nsim + 1), including
ties, to avoid exact zero and one. The maintainer and contact are unchanged.

## Executed checks

- R 4.6.1, x86_64 macOS Sonoma 14.8.7, current indexed mvtnorm 1.4-2,
  numDeriv 2016.8-1.1, and testthat 3.3.2.
- Build succeeded. All 17 testthat cases and 174 expectations passed without
  warnings or skips. Namespace isolation and examples passed.
- Supplementary R CMD check --as-cran completed with Status: OK, including
  PDF/HTML manuals, with incoming feasibility and remote clock checks disabled.
- Normal incoming feasibility stopped because R could not resolve repository
  hosts. This limitation is not treated as a successful full check.
- The fresh CRAN package index listed no direct reverse dependencies on fgof
  in Depends, Imports, LinkingTo, Suggests, or Enhances.
- Paired numerical comparisons covered 75 continuous-model datasets and twelve
  correlated normal/t datasets, using analytical derivatives and identical
  multiplier draws. These are numerical validation pilots, not claims of exact
  finite-sample calibration.

Evidence and complete selected check logs are under
misc/audit/derivative-integration/. Methodology and limits are described in
analysis/derivative-scale/INTEGRATION.md; development files are excluded from
source package builds.

## Before submission

Complete normal online incoming checks and cross-platform CI on the final
commits, refresh these results if the source changes, and obtain the maintainer's
explicit approval to submit. No merge, release, or submission has occurred.

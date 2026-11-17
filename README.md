# fgof

`fgof` tests goodness of fit for continuous distributions with estimated
parameters using multiplier or parametric bootstrap procedures. This repository
revives CRAN release 0.2-1, archived on 2022-04-27.

```r
library(fgof)
set.seed(123)
x <- rlnorm(100, meanlog = 1, sdlog = 0.5)
fit <- gof.test(x, "lnorm", nsim = 199,
                start = list(meanlog = 1, sdlog = 0.5))
fit$p.value
```

The default uses centered normal multipliers. Set `simulation = "pb"` for
parametric bootstrap with refitting. See `help("gof.test")` for statistics,
parameter constraints, derivative callbacks, and optimizer controls. The
starting values specify parameters to estimate, rather than known parameters.
Bootstrap p-values use `(count + 0.5) / (nsim + 1)`, where `count` includes ties,
so they remain strictly between zero and one.

## Development

Install the R dependencies `mvtnorm` and `numDeriv`, plus `testthat` for tests,
then run:

```sh
make build    # Create fgof_<version>.tar.gz
make check    # Build and run R CMD check --as-cran
make install  # Build and install in the selected R library
```

Regression tests use `testthat` and run as part of `R CMD check`. To run them
during development, use `Rscript --vanilla -e 'testthat::test_local()'`.
A separate base R check verifies namespace isolation. The GitHub Actions
workflow checks macOS, Windows, and Linux, including R-devel and oldrel-1 after
the repository is pushed. No remote check has been run from this workspace.

## Revival records

- [Initial source audit](REVIVAL_AUDIT.md)
- [Revival progress and remaining requirements](REVIVAL_REPORT.md)
- [Recovery provenance and original source](misc/recovery/README.md)
- [Methodological decisions](misc/validation/methodology.md)
- [Validation scripts and results](misc/validation/README.md)
- [Changes](NEWS.md)

The original statistical procedures are described by
[Kojadinovic and Yan (2012)](https://doi.org/10.1002/cjs.11135).
The package is licensed under GPL (>= 3).

## Numerical differentiation

With `method.args` omitted, the multiplier bootstrap checks derivatives for
standard and user-defined models using the same generic procedure. It infers
local parameter scales from log-density changes or uses an existing optimizer
`control$parscale`. Steps are limited by a relative cap and tested feasible
radii, then checked for stability separately for the score and CDF. Functions
must give deterministic values. Unreliable derivatives produce an error with
suggested alternatives.

Explicit `method.args` retains the previous numDeriv behavior and bypasses the
automatic checks. For example, `list(d = 0, eps = unname(h), zero.tol = Inf,
r = 6)` supplies absolute initial steps. Give exactly one positive entry in `h`
per estimated parameter in density argument order; names do not align entries.
Check sensitivity to plausible steps, since smaller steps can lose precision.
Analytical callbacks remain available and must incorporate fixed parameters.
These changes leave the score-covariance and influence-matrix formulas intact.

See [the derivative integration report](analysis/derivative-scale/INTEGRATION.md)
for executed comparisons, limitations, and package-check results.

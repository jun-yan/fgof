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
workflow checks macOS, Windows, and Linux, including R-devel and oldrel-1 when
changes are pushed.

## Repository scope

This repository contains the R package, documentation, routine tests, and CI.
Supporting investigations, benchmarks, recovery archives, and manuscript assets
are maintained separately under the untracked local `investigation/` directory.
See [NEWS.md](NEWS.md) for package changes.

The original statistical procedures are described by
[Kojadinovic and Yan (2012)](https://doi.org/10.1002/cjs.11135).
The package is licensed under GPL (>= 3).

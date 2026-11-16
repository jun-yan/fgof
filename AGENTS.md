# Repository Guidelines

## Project Structure & Module Organization

`fgof` is an R package being revived from CRAN release 0.2-1. `R/gof.R`
contains `gof.test()` and fitting/simulation helpers; `R/derMvt.R` contains
multivariate t derivatives; `R/mvdwrapper.R` provides distribution wrappers.
`man/gof.Rd` documents the public API. `DESCRIPTION` declares metadata and
dependencies; `NAMESPACE` currently exports only `gof.test`.

`misc/recovery/` holds the original tarball, historical CRAN checks, provenance,
and a Git bundle. Preserve these artifacts. No test suite, Makefile, CI workflow,
or pkgdown configuration exists yet.

## Build, Test, and Development Commands

Run commands from the repository root:

- `Rscript --vanilla -e 'install.packages(c("mvtnorm", "numDeriv"), repos="https://cloud.r-project.org")'`: install declared dependencies.
- `R CMD build .`: build a source tarball.
- `R CMD check --as-cran fgof_0.2-1.tar.gz`: check the built package; adjust the filename after version changes.
- `R CMD INSTALL fgof_0.2-1.tar.gz`: install the built package locally.

Before routine builds, add `.Rbuildignore` entries for `misc` and this guide.
Historical failures remain unresolved; report check results accurately.

## Coding Style & Naming Conventions

Use two-space indentation, `<-` for assignment, and spaces around operators.
Follow surrounding code and preserve existing function and argument names.
Keep edits focused; avoid unrelated reformatting. Update Rd documentation when
changing public behavior. No formatter, linter, or roxygen workflow is configured.
Declare dependency functions explicitly through namespace imports or `pkg::fun`.

## Testing Guidelines

No testing framework or coverage threshold is configured. Validate changes with
`R CMD check --as-cran`, including examples. Add focused regression tests for
repairs; use `tests/*.R` with `stopifnot()` until a framework is adopted.
Use descriptive filenames such as `tests/qgrid.R`. Set seeds for simulations,
use appropriate numerical tolerances, and keep simulation counts practical.
Prioritize the archived `qgrid` example failure and missing namespace imports.

## Commit & Pull Request Guidelines

The preserved history contains one commit, `version 0.2-1`, so no broader message
convention is established. Use concise imperative messages describing one change.
Keep recovery artifacts separate from repairs. PRs should explain the problem,
behavior changes, validation commands/results, and any remaining failures. Link
relevant issues and identify statistical or API compatibility implications.

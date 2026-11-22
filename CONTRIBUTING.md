# Contributing to fgof

## Build and test

Run `make build` from the repository root to create the source tarball and
`make check` to run `R CMD check --as-cran`, including examples, tests, and
vignettes. The [README](README.md) describes the development workflow.

Add focused regression tests under `tests/testthat/` when fixing bugs. Set
seeds for simulations, use appropriate numerical tolerances, and keep
simulation counts practical.

## Code and documentation

Use two-space indentation, `<-` for assignment, and spaces around operators.
Follow the surrounding code and keep edits focused. Declare dependency
functions explicitly through namespace imports or `pkg::fun`. When changing
public behavior, update the roxygen documentation and regenerate the help
files with `make document`.

## Commits and pull requests

Use concise, imperative commit messages. In pull requests, explain the problem,
the behavior that changed, and the checks you ran. Link relevant issues and
describe any statistical or API compatibility implications.

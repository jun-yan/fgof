# fgof development version

## Compatibility

- Normalize `simulation` and `qgrid` choices before scalar conditions, restoring
  their documented defaults on current R. Invalid choices now report an error.
- Import functions used from stats, mvtnorm, and numDeriv so namespace-only calls
  work without attaching dependency packages.
- Exclude recovery records and development files from source package builds.
- Remove the recovered tarball's generated MD5 manifest from the maintained
  source tree; package builds generate a fresh manifest.
- Add dependency-free, fixed-seed regression tests.

# Changelog

All notable changes to the `dcebe` Python package are documented here.
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
and the [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) layout.

The MATLAB reference implementation is tracked separately in this same
repo's git history; algorithmic semantics are unchanged across the port.

## [Unreleased]

## [1.0.0] — 2026-05-11

First stable release.

The MATLAB reference DCEBE has been published since 2018 (paper:
Bendinger, Debus, Glowa, Karger, Peter, Storath, *Bolus arrival time
estimation in dynamic contrast-enhanced MRI of small animals based on
spline models*, Physics in Medicine & Biology 64(4), 2019). The Python
re-implementation in this repository was developed in 2026 as a port
of the MATLAB code; the user-facing API (`estimate_bat`) is unchanged
from the 0.1.0 dev release.

### Added

- Python re-implementation of `DCEBE_estimateBAT.m` as `estimate_bat`,
  returning an `EstimateResult` dataclass with 1-based MATLAB-compatible
  sample units.
- `pip install dcebe` from PyPI. Pure Python (NumPy + SciPy); no
  compiled wheels, no Rust, no MATLAB runtime required.
- Auto-create GitHub Release on tag push (alongside PyPI publish).
- `CITATION.cff` with `version` and `date-released` fields.
- README harmonised with the lab-repo family (badge block,
  Quickstart-first ordering, "See also" section linking the five
  sibling repos, License footer).

### Changed

- `LICENSE` cleaned up: non-standard "if this code is used in any
  scientific publication, the following paper shall be cited" clause
  removed. The LICENSE file is now pure Apache-2.0 matching the
  declaration in `CITATION.cff` and `pyproject.toml`. The citation
  request lives in `CITATION.cff` and the README.

[Unreleased]: https://github.com/mstorath/DCEBE/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/mstorath/DCEBE/releases/tag/v1.0.0

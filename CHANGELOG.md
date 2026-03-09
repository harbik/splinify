# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.1] - 2026-03-07

### Fixed

- `serde_json` serialization of `SplineCurve` in tests and examples now correctly wraps the result in `SplineCurveData`.
- README code examples marked `no_run` so they are not executed as doc tests (they require external CSV files).
- Corrected type name in LED spectral data README example (`SplineCurveFit::<3>` → `CubicSplineFit`).
- Fixed `repository` URL in `Cargo.toml` (was pointing to `dierckx-sys`, now points to `splinify`).

## [0.2.0] - 2025-12-14

### Changed

- `plotters` is now an **optional** dependency, enabled via the `plot` feature (`features = ["plot"]`).
  Previously it was always compiled in, requiring its full dependency tree for all users.
- `spliny` and `dierckx-sys` dependencies updated from git path references to published crates.io versions.

## [0.1.0] - 2021-12-27

Initial release on [crates.io](https://crates.io/crates/splinify).

### Added

- `SplineCurveFit<K>` — wraps Dierckx' `curfit` subroutine to fit a degree-*K* spline to *(x, y)* data.
- `ParameterSplineCurveFit<K, N>` — wraps Dierckx' `concur` subroutine for optionally constrained,
  parametric N-dimensional curve fitting.
- Fit methods on both types:
  - `smoothing_spline(rms)` — automatically selects knots to achieve a target RMS error.
  - `smoothing_spline_optimize` — iterative optimizer variant for `ParameterSplineCurveFit`.
  - `cardinal_spline(dt)` — weighted least-squares spline with equidistant knots spaced `dt` apart.
  - `interpolating_spline()` — exact-fit spline with knots at the data locations.
- `begin_constraints` / `end_constraints` on `ParameterSplineCurveFit` for specifying derivative
  boundary conditions at the endpoints.
- `SplineCurveData` — serializable wrapper around `SplineCurve` (implements `serde::Serialize`).
- Type aliases for common degree/dimension combinations:
  `LinearSplineFit`, `CubicSplineFit`, `QuinticSplineFit`,
  `CubicSplineFit1D` / `2D` / `3D`, `QuinticSplineFit1D` / `2D` / `3D`, etc.
- CSV helpers `read_csv_xy` and `read_csv_uxy`.
- Examples: `lissajous`, `bb-locus`, `pezzack`, `led-fit`.
- Apache-2.0 / MIT dual license.

[0.2.1]: https://github.com/harbik/splinify/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/harbik/splinify/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/harbik/splinify/releases/tag/v0.1.0

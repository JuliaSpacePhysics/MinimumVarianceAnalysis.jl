# Changelog

## [Unreleased]

### Changed

- **Breaking**: `mva_eigen` and `mva` find the time dimension of containers through TimeseriesUtilities (`DimArray`, `DataVariable`); plain matrices without `dim` throw `ArgumentError` instead of assuming samples along rows.
- **Breaking**: `mva(V, F::Eigen)` replaces `transform(V, F)` and its `rotate` alias; `mva` output has anonymous component axes, since components then lie along the eigenvectors.
- More than 3 components use a dynamic variance matrix instead of a static one.
- **Breaking**: `mva_eigen` returns `MVAResult` instead of `Eigen`. `normal(F)` and the quality check take `field` from `F`.
- **Breaking**: `check_mva_eigen(F; field, verbose)` is replaced by `issuccess(F; r0)`, a pure predicate; `mva_eigen(...; check = true)` warns when it fails.
- **Breaking**: `B_x3_error` and `E_x1_error` are replaced by `normal_component_error(F)`.
- **Breaking**: `mva_eigen` drops the `sort` keyword; eigenvalues are always λ₁ ≥ λ₂ ≥ λ₃.

# Changelog

## [Unreleased]

### Changed

- **Breaking**: `mva_eigen`, `mva` and `convection_efield` find the time dimension of containers through TimeseriesUtilities (`DimArray`, `DataVariable`); plain matrices without `dim` throw `ArgumentError` instead of assuming samples along rows.
- **Breaking**: `mva(V, F::Eigen)` replaces `transform(V, F)` and its `rotate` alias; `mva` output has anonymous component axes, since components then lie along the eigenvectors.
- More than 3 components use a dynamic variance matrix instead of a static one.

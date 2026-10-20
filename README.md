# MinimumVarianceAnalysis

The main purpose of minimum or maximum variance analysis (MVA) is to find an estimator for the direction normal $\hat{𝐧}$ to one-dimensional structure, by minimisation of

$$
σ^2=\frac{1}{M} ∑_{m=1}^M | (𝐁^{(m)}-⟨𝐁⟩) · \hat{𝐧}|^2.
$$

## Quickstart

```julia
using Pkg; Pkg.add("MinimumVarianceAnalysis")
using MinimumVarianceAnalysis
using MinimumVarianceAnalysis: normal, normal_component_error

Bl = cosd.(0:30:180)
Bm = sind.(0:30:180)
Bn = 0.1 .* ones(7)
B = hcat(Bl, Bm, Bn)  # rows are samples

F = mva_eigen(B; dim = 1)    # MVAResult: values (λ₁ ≥ λ₂ ≥ λ₃), eigenvector columns `vectors`, mean, n, field
normal(F)                    # minimum-variance direction for magnetic-field MVA
issuccess(F)                 # normal well determined if |λ₂/λ₃| > r0 (default 5)
normal_component_error(F)    # statistical error of ⟨B⟩·normal
mva(B; dim = 1)              # rotate samples into variance frame
mva(B, F; dim = 1)           # rotate any data (e.g. velocity) into F's frame
```

For electric-field maximum variance analysis, pass `field = :E`:

```julia
E = hcat(0.1 .* ones(7), 0.05 .* ones(7), cosd.(0:30:180))
F = mva_eigen(E; dim = 1, field = :E)
normal(F)    # maximum-variance direction; issuccess uses |λ₁/λ₂|
```

Containers that carry their times (`DimArray`, SpaceDataModel's `DataVariable`) default to their time dimension.
Unitful, DynamicQuantities and FlexUnits inputs are supported:

```julia
using Unitful

F = mva_eigen(B .* u"nT"; dim = 1)
F.values    # eigenvalues carry squared field units
F.vectors   # eigenvectors are dimensionless directions
```

See [SPEDAS](https://juliaspacephysics.github.io/SPEDAS.jl/dev/explanations/coords/) for more details. See [SPEDAS validation](https://juliaspacephysics.github.io/SPEDAS.jl/dev/validation/pyspedas/#Minimum-variance-analysis) for comparison with a Python implementation (pyspedas).

**Documentation**: [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSpacePhysics.github.io/MinimumVarianceAnalysis.jl/dev/)

## Reference

- [Sonnerup, B. U. Ö., & Scheible, M. (1998). Minimum and maximum variance analysis. ISSI Scientific Reports Series, 1, 185–220.](https://ui.adsabs.harvard.edu/abs/1998ISSIR...1..185S/abstract)

## Features and Roadmap

- [x] Maximum Variance Analysis on Magnetic Field (MVAB)
    - [ ] Constraint $〈B₃〉 = 0$
- [x] Maximum Variance Analysis on Electric Field (MVAE)
    - [ ] Transformation to/from moving frame of reference
- [ ] Minimum Variance Analysis on Mass Flux (MVAρv)
- [ ] Application to 2-D Structures

## Notes

- Anisotropic fluctuations have been shown to lead to larger errors in normal determinations.
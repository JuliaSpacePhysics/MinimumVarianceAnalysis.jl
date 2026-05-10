# MinimumVarianceAnalysis

The main purpose of minimum or maximum variance analysis (MVA) is to find an estimator for the direction normal $\hat{𝐧}$ to one-dimensional structure, by minimisation of

$$
σ^2=\frac{1}{M} ∑_{m=1}^M | (𝐁^{(m)}-⟨𝐁⟩) · \hat{𝐧}|^2.
$$

## Quickstart

```julia
using Pkg; Pkg.add("MinimumVarianceAnalysis")
using MinimumVarianceAnalysis
using MinimumVarianceAnalysis: normal

Bl = cosd.(0:30:180)
Bm = sind.(0:30:180)
Bn = 0.1 .* ones(7)
B = hcat(Bl, Bm, Bn)

F = mva_eigen(B)
normal(F)          # minimum-variance direction for magnetic-field MVA
mva(B)             # rotate samples into variance frame
```

For electric-field maximum variance analysis, pass `field = :E`:

```julia
E = hcat(0.1 .* ones(7), 0.05 .* ones(7), cosd.(0:30:180))
F = mva_eigen(E; field = :E)
normal(F; field = :E)
```

Unitful, DynamicQuantities, FlexUnits, and DimensionalData inputs are supported:

```julia
using Unitful

F = mva_eigen(B .* u"nT")
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

## Status

[![DOI](https://zenodo.org/badge/1046063112.svg)](https://doi.org/10.5281/zenodo.18635364)
[![version](https://juliahub.com/docs/General/MinimumVarianceAnalysis/stable/version.svg)](https://juliahub.com/ui/Packages/General/MinimumVarianceAnalysis)

[![Build Status](https://github.com/JuliaSpacePhysics/MinimumVarianceAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSpacePhysics/MinimumVarianceAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/MinimumVarianceAnalysis.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/MinimumVarianceAnalysis.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

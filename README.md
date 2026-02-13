# MinimumVarianceAnalysis

[![DOI](https://zenodo.org/badge/1046063112.svg)](https://doi.org/10.5281/zenodo.18635364)
[![version](https://juliahub.com/docs/General/MinimumVarianceAnalysis/stable/version.svg)](https://juliahub.com/ui/Packages/General/MinimumVarianceAnalysis)


[![Build Status](https://github.com/JuliaSpacePhysics/MinimumVarianceAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSpacePhysics/MinimumVarianceAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/MinimumVarianceAnalysis.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/MinimumVarianceAnalysis.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

The main purpose of minimum or maximum variance analysis (MVA) is to find an estimator for the direction normal $\hat{𝐧}$ to an approximately one-dimensional structure, by minimisation of

$$
σ^2=\frac{1}{M} ∑_{m=1}^M | (𝐁^{(m)}-⟨𝐁⟩) · \hat{𝐧}|^2.
$$

See [SPEDAS](https://juliaspacephysics.github.io/SPEDAS.jl/dev/explanations/coords/) for more details. See [SPEDAS validation](https://juliaspacephysics.github.io/SPEDAS.jl/dev/validation/pyspedas/#Minimum-variance-analysis) for comparison with a Python implementation (pyspedas).

**Installation**: at the Julia REPL, run `using Pkg; Pkg.add("MinimumVarianceAnalysis")`

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
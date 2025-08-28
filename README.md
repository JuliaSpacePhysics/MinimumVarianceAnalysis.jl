# MinimumVarianceAnalysis

[![Build Status](https://github.com/JuliaSpacePhysics/MinimumVarianceAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSpacePhysics/MinimumVarianceAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/MinimumVarianceAnalysis.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/MinimumVarianceAnalysis.jl)

The main purpose of minimum or maximum variance analysis (MVA) is to find an estimator for the direction normal $\hat{𝐧}$ to an approximately one-dimensional structure, by minimisation of

$$
σ^2=\frac{1}{M} ∑_{m=1}^M | (𝐁^{(m)}-⟨𝐁⟩) · \hat{𝐧}|^2.
$$

See [SPEDAS](https://juliaspacephysics.github.io/SPEDAS.jl/dev/explanations/coords/) for more details. See [SPEDAS validation](https://juliaspacephysics.github.io/SPEDAS.jl/dev/validation/pyspedas/#Minimum-variance-analysis) for comparison with a Python implementation (pyspedas).

## Reference

- [Sonnerup, B. U. Ö., & Scheible, M. (1998). Minimum and maximum variance analysis. ISSI Scientific Reports Series, 1, 185–220.](https://ui.adsabs.harvard.edu/abs/1998ISSIR...1..185S/abstract)

## Roadmap

- [ ] Minimum Variance Analysis on Mass Flux (MVAρv)
- [ ] Maximum Variance Analysis on Electric Field (MVAE)
- [ ] Application to 2-D Structures
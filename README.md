# MinimumVarianceAnalysis

The main purpose of minimum or maximum variance analysis (MVA) is to find an estimator for the direction normal $\hat{𝐧}$ to one-dimensional structure, by minimisation of the variance along $\hat{𝐧}$ in excess of background noise

$$
σ^2 = \frac{1}{M} ∑_{m=1}^{M} |(𝐁^{(m)}-⟨𝐁⟩) · \hat{𝐧}|^2 - \frac{1}{∑_k M_k} ∑_k ∑_{m=1}^{M_k} |(𝐁_k^{(m)}-⟨𝐁_k⟩) · \hat{𝐧}|^2,
$$

where $𝐁_k$ are optional quiet intervals of $M_k$ samples.

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
issuccess(F; r0 = 5)         # normal well determined if excess variance along the intermediate axis > (r0 - 1) × variance along the normal
normal_component_error(F; M = F.n)  # statistical error of field·normal
mva(B; dim = 1)              # rotate samples into variance frame
mva(B, F; dim = 1)           # rotate any data (e.g. velocity) into F's frame
```

Pass the quiet intervals $𝐁_k$ as `background`:

```julia
F = mva_eigen(B[transition, :]; dim = 1, background = (B[before, :], B[after, :]))
F.values     # excess variance over the background (may be negative)
F.noise      # background covariance in the eigenbasis
```

Anisotropic noise (e.g. Alfvénic fluctuations, with little variance along 𝐁₀) biases the normal toward the noise's least-variance direction. The bias is systematic, so `issuccess` and `normal_component_error`, which measure only random error, can still report a reliable normal ([test](test/runtests.jl)). Pass `background` to remove the bias.

- Default to a background whenever quiet data ≳ 10× longer than the transition is available; under isotropic noise it then costs little.
- A shorter background adds noise: use it only if its own eigenvalue ratio shows clear anisotropy (≳ 2), or if normals with and without it disagree beyond `normal_component_error`.
- Take intervals on both sides, outside the transition, in the same plasma regime, free of other structures and trends.

For electric-field maximum variance analysis, pass `field = :E`:

```julia
E = hcat(0.1 .* ones(7), 0.05 .* ones(7), cosd.(0:30:180))
F = mva_eigen(E; dim = 1, field = :E)
normal(F)    # maximum-variance direction
```

Containers that carry their times (`DimArray`, SpaceDataModel's `DataVariable`) default to their time dimension.
Unitful, DynamicQuantities and FlexUnits inputs are supported:

```julia
using Unitful

F = mva_eigen(B .* u"nT"; dim = 1)
F.values    # eigenvalues carry squared field units
F.vectors   # eigenvectors are dimensionless directions
```

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

- Variance outside the transition is noise: MVA over the transition itself beat the customary ±1.5 transition widths. Locate transitions with [ScaleSpaceFeatures.jl](https://github.com/JuliaSpacePhysics/ScaleSpaceFeatures.jl).
- Correlated samples (turbulence, filtering) might make `F.n` overstate the information; pass independent samples `M ≈ 2 f_c T` (cutoff `f_c`, window `T`) to `normal_component_error`.
- If 𝐁 varies along one direction (λ₂ ≈ λ₃), one spacecraft only constrains 𝐧 ⟂ Δ𝐁; no variance variant recovers it. Use multi-spacecraft gradients or timing ([MultiSpacecraftAnalysis.jl](https://github.com/JuliaSpacePhysics/MultiSpacecraftAnalysis.jl)).

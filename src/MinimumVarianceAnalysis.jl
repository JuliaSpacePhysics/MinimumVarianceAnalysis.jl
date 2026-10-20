# References:
# - https://github.com/henry2004y/VisAnaJulia/blob/master/src/MVA.jl
# - https://pyspedas.readthedocs.io/en/latest/coords.html#pyspedas.minvar
"""
Estimate the direction normal ``𝐧̂`` to a one-dimensional structure using minimum or maximum variance analysis (MVA).

Main Functions: [`mva`](@ref), [`mva_eigen`](@ref), [`normal`](@ref)
"""
module MinimumVarianceAnalysis

using LinearAlgebra
using StaticArrays
using TimeseriesUtilities: resolve_dim, unwrap, rebuild, mapsamples
export mva, mva_eigen, check_mva_eigen
export convection_efield

const SV3 = SVector{3}

include("utils.jl")

"""
    mva(V, F=V; dim=nothing, kwargs...)
    mva(V, F::Eigen; dim=nothing)

Transform time series `V` into the variance frame (LMN) of reference field `F`, or of an [`mva_eigen`](@ref) result `F`.
`kwargs` go to [`mva_eigen`](@ref).
"""
mva(V, F = V; dim = nothing, kwargs...) = mva(V, mva_eigen(F; dim, kwargs...); dim)

function mva(V, F::Eigen; dim = nothing)
    d = resolve_dim(V, dim)
    data = unwrap(V)
    out = d == 1 ? data * F.vectors : F.vectors' * data
    # components now lie along the eigenvectors, so the input's component labels no longer apply
    c = 3 - d
    return rebuild(V, out, Base.setindex((Val(1), Val(2)), nothing => axes(out, c), c))
end

@inline function sorteigen(F; sortby = abs, rev = true)
    order = sortperm(F.values; rev, by = sortby)
    return Eigen(F.values[order], F.vectors[:, order])
end

# Static path for the 3-component case only: a literal size leaves `M` a small Union the compiler splits, so no dynamic dispatch.
function _mva_eigen(B; sort = (;))
    M = size(B, 2) == 3 ? _variance_matrix(B, Val(3)) : _variance_matrix(B)
    return sorteigen(eigen(M); sort...)
end

function _variance_matrix(B, v::Val{N}) where {N}
    n = size(B, 1)
    B̄ = SVector{N}(sum(Bc) / n for Bc in eachcol(B))
    # Center before the products: the raw-moment form ⟨BBᵀ⟩ - B̄B̄ᵀ cancels catastrophically when |B̄| ≫ δB.
    M = zero(SMatrix{N, N, eltype(B̄)})
    @inbounds for k in axes(B, 1)
        d = SVector(ntuple(i -> B[k, i], v)) - B̄
        M += d * d'
    end
    return Symmetric(M / n)
end

function _variance_matrix(B)
    C = B .- sum(B; dims = 1) ./ size(B, 1)
    return Symmetric(C' * C / size(B, 1))
end

"""
    mva_eigen(x::AbstractMatrix; dim = nothing, sort=(;), check=false) -> F::Eigen

Perform minimum variance analysis of the magnetic field `B` or maximum variance analysis of the electric field `E` when `field=:E`.

Samples lie along `dim`, by default the time dimension of a container (`DimArray`, `DataVariable`, …); plain matrices must pass `dim`.

Return `Eigen` factorization object `F` which contains the eigenvalues in `F.values` and the eigenvectors in the columns of the matrix `F.vectors`. The `k`th eigenvector can be obtained from the slice `F.vectors[:, k]`.

Set `check=true` to check the reliability of the result.

## Notes

For a one-dimensional current layer, the tangential electric field components are approximately constant across the boundary, while the normal component exhibits the largest variation. Therefore, the eigenvector corresponding to the **maximum eigenvalue** ``λ_1`` (first column of `F.vectors`) gives an estimate of the boundary normal direction.
"""
function mva_eigen(B; dim = nothing, sort = (;), check = false, field = :B)
    data = unwrap(B)
    F = _mva_eigen(resolve_dim(B, dim) == 1 ? data : transpose(data); sort)
    check && check_mva_eigen(F; field)
    return F
end

function _field(dimension, electric_dimension, magnetic_dimension)
    dimension == electric_dimension && return :E
    dimension == magnetic_dimension && return :B
    return nothing
end

function _quantity_mva_eigen(B, scale, raw_values, infer_field_dimension; field = nothing, kwargs...)
    field = @something field infer_field_dimension(scale) :B
    F = mva_eigen(raw_values(B); field, kwargs...)
    return Eigen(F.values .* scale^2, F.vectors)
end

"""
    check_mva_eigen(F; r0=5, verbose=false, field = :B)

Check the quality of the MVA result.

If λ₁ ≥ λ₂ ≥ λ₃ are 3 eigenvalues of the constructed matrix M. For MVAB, a good indicator of nice results should have ``|λ₂ / λ₃| > r₀`` (default ``r₀ = 5``).

For MVAE, a reliable normal direction requires the maximum eigenvalue ``λ₁`` to be well-separated from the intermediate eigenvalue ``λ₂``. The ratio ``|λ₁ / λ₂| > r₀`` is used as a quality indicator.
"""
function check_mva_eigen(F; r0 = 5, verbose = false, field = :B)
    r = field == :E ? abs(F.values[1] / F.values[2]) : abs(F.values[2] / F.values[3])
    flag = r > r0
    verbose && begin
        println(F.vectors)
        if field == :E
            println("Ratio of maximum variance to intermediate variance = ", r)
        else
            println("Ratio of intermediate variance to minimum variance = ", r)
        end
        flag && @info "Seems to be a proper MVA attempt!"
    end
    flag || @warn "Take the MVA result with a grain of salt!"
    return flag
end

################
# Error Estimate
################

"""
    Δφij(λᵢ, λⱼ, λ₃, M)

Calculate the phase error between components i and j according to:
|Δφᵢⱼ| = |Δφⱼᵢ| = √(λ₃/(M-1) * (λᵢ + λⱼ - λ₃)/(λᵢ - λⱼ)²)

Parameters:

  - λᵢ: eigenvalue i
  - λⱼ: eigenvalue j
  - λ₃: smallest eigenvalue (λ₃)
  - M: number of samples
"""
function Δφij(λᵢ, λⱼ, λ₃, M)
    return sqrt((λ₃ / (M - 1)) * (λᵢ + λⱼ - λ₃) / (λᵢ - λⱼ)^2)
end

include("mvab.jl")
include("mvae.jl")
include("workload.jl")

end

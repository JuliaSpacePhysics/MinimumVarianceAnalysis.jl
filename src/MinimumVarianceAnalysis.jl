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
export mva, mva_eigen, check_mva_eigen
export convection_efield

const SV3 = SVector{3}

include("utils.jl")
include("unitful.jl")

"""
    mva(V, F=V; dim=nothing, kwargs...)

Transform a timeseries `V` into the LMN coordinate system based on the minimum/maximum variance analysis of reference field `F` along the `dim` dimension (time).

See also: [`mva_eigen`](@ref), [`transform`](@ref)
"""
function mva(V, F = V; dim = nothing, kwargs...)
    E = mva_eigen(F; dim, kwargs...)
    return transform(V, E; dim)
end

@inline function sorteigen(F; sortby = abs, rev = true)
    order = sortperm(F.values; rev, by = sortby)
    return Eigen(F.values[order], F.vectors[:, order])
end

function _mva_eigen(B, ::Val{N}; sort = (;)) where {N}
    n = size(B, 1)
    Bslices = eachcol(B)
    B̄ = SVector{N}(sum(Bc) / n for Bc in Bslices)
    M = SMatrix{N, N}(dot(Bi, Bj) / n for Bi in Bslices, Bj in Bslices) - B̄ * B̄'
    return sorteigen(eigen(M); sort...)
end

"""
    mva_eigen(x::AbstractMatrix; dim = nothing, sort=(;), check=false) -> F::Eigen

Perform minimum variance analysis of the magnetic field `B` or maximum variance analysis of the electric field `E` when `field=:E`.

`x` varies along the `dim` dimension.

Return `Eigen` factorization object `F` which contains the eigenvalues in `F.values` and the eigenvectors in the columns of the matrix `F.vectors`. The `k`th eigenvector can be obtained from the slice `F.vectors[:, k]`.

Set `check=true` to check the reliability of the result.

## Notes

For a one-dimensional current layer, the tangential electric field components are approximately constant across the boundary, while the normal component exhibits the largest variation. Therefore, the eigenvector corresponding to the **maximum eigenvalue** ``λ_1`` (first column of `F.vectors`) gives an estimate of the boundary normal direction.
"""
function mva_eigen(B; dim = nothing, sort = (;), check = false, field = :B)
    dim = something(dim, 1)
    in = dim == 1 ? B : B'
    N = size(in, 2)
    F = _mva_eigen(in, Val(N); sort)
    check && check_mva_eigen(F; field)
    return F
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

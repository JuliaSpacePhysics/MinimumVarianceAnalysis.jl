# References:
# - https://github.com/henry2004y/VisAnaJulia/blob/master/src/MVA.jl
# - https://pyspedas.readthedocs.io/en/latest/coords.html#pyspedas.minvar
module MinimumVarianceAnalysis

using LinearAlgebra
using StaticArrays
using Unitful: Quantity, ustrip, unit
export mva, mva_eigen, check_mva_eigen

const SV3 = SVector{3}

include("utils.jl")

"""
    mva(V, B=V; dim=nothing, kwargs...)

Rotate a timeseries `V` into the LMN coordinates based on the reference field `B` along the `dim` dimension (time).

See also: [`mva_eigen`](@ref), [`rotate`](@ref)
"""
function mva(V, B = V; dim = nothing, kwargs...)
    E = mva_eigen(B; dim, kwargs...)
    return rotate(V, E; dim)
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
    mva_eigen(B::AbstractMatrix; dim = nothing, sort=(;), check=false) -> F::Eigen

Perform minimum variance analysis for `B` (which varies along the `dim` dimension).
    
Return `Eigen` factorization object `F` which contains the eigenvalues in `F.values` and the eigenvectors in the columns of the matrix `F.vectors`.

Set `check=true` to check the reliability of the result.

The `k`th eigenvector can be obtained from the slice `F.vectors[:, k]`.
"""
function mva_eigen(B; dim = nothing, sort = (;), check = false)
    dim = something(dim, 1)
    in = dim == 1 ? B : B'
    N = size(in, 2)
    F = _mva_eigen(in, Val(N); sort)
    check && check_mva_eigen(F)
    return F
end

function mva_eigen(B::AbstractMatrix{Q}; kwargs...) where {Q <: Quantity}
    F = mva_eigen(ustrip(B); kwargs...)
    return Eigen(F.values * unit(Q)^2, F.vectors)
end


"""
    check_mva_eigen(F; r=5, verbose=false)

Check the quality of the MVA result.

If λ₁ ≥ λ₂ ≥ λ₃ are 3 eigenvalues of the constructed matrix M, then a good
indicator of nice fitting LMN coordinate system should have ``abs(λ₂ / λ₃) > r``.
"""
function check_mva_eigen(F; r0 = 5, verbose = false)
    r = abs(F.values[2] / F.values[3])
    flag = r > r0
    verbose && begin
        println(F.vectors)
        println("Ratio of intermediate variance to minimum variance = ", r)
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

"""
Calculate the composite statistical error estimate for ⟨B·x₃⟩:
|Δ⟨B·x₃⟩| = √(λ₃/(M-1) + (Δφ₃₂⟨B⟩·x₂)² + (Δφ₃₁⟨B⟩·x₁)²)

Parameters:

  - λ₁, λ₂, λ₃: eigenvalues in descending order
  - M: number of samples
  - B: mean magnetic field vector
  - x₁, x₂, x₃: eigenvectors
"""
function B_x3_error(λ₁, λ₂, λ₃, M, B, x₁, x₂, x₃)
    Δφ₃₂ = Δφij(λ₃, λ₂, λ₃, M)
    Δφ₃₁ = Δφij(λ₃, λ₁, λ₃, M)
    B_x₂ = dot(B, x₂)
    B_x₁ = dot(B, x₁)
    return sqrt(λ₃ / (M - 1) + (Δφ₃₂ * B_x₂)^2 + (Δφ₃₁ * B_x₁)^2)
end

B_x3_error(F::Eigen, M, B) = B_x3_error(F.values..., M, B, eachcol(F.vectors)...)

include("workload.jl")

end

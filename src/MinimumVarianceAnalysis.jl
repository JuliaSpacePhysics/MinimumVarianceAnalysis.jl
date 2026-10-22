module MinimumVarianceAnalysis

using LinearAlgebra
import LinearAlgebra: issuccess
using StaticArrays
using TimeseriesUtilities: resolve_dim, unwrap, rebuild
export mva, mva_eigen, issuccess

struct MVAResult{V, M, B}
    values::V
    vectors::M
    mean::B
    n::Int
    field::Symbol
end

# The error and quality formulas assume a 3-D field: λ₃ is the minimum eigenvalue.
_check3(F::MVAResult) = length(F.values) == 3 || throw(ArgumentError("defined for 3 components, got $(length(F.values))"))

function Base.show(io::IO, mime::MIME"text/plain", F::MVAResult)
    println(io, "MVAResult (field = :", F.field, ", n = ", F.n, ")")
    println(io, "values:")
    show(io, mime, F.values)
    println(io, "\nvectors:")
    return show(io, mime, F.vectors)
end

@views normal(F::MVAResult) = F.field == :E ? F.vectors[:, 1] : F.vectors[:, end]

is_right_handed(F::MVAResult) = det(F.vectors) > 0

"""
    mva(V, F=V; dim=nothing, kwargs...)
    mva(V, F::MVAResult; dim=nothing)

Transform time series `V` into the variance frame (LMN) of reference field `F`, or of an [`mva_eigen`](@ref) result `F`.
"""
mva(V, F = V; dim = nothing, kwargs...) = mva(V, mva_eigen(F; dim, kwargs...); dim)

function mva(V, F::MVAResult; dim = nothing)
    d = resolve_dim(V, dim)
    data = unwrap(V)
    out = d == 1 ? data * F.vectors : F.vectors' * data
    # components now lie along the eigenvectors, so the input's component labels no longer apply
    c = 3 - d
    return rebuild(V, out, Base.setindex((Val(1), Val(2)), nothing => nothing, c))
end

# Static path for the 3-component case only
function _mva_eigen(B, field)
    return size(B, 2) == 3 ? _mva_eigen(B, Val(3), field) : _mva_eigen(B, nothing, field)
end

function _mva_eigen(B, N, field)
    B̄, M = _moments(B, N)
    λ, V = eigen(M)
    order = sortperm(λ; rev = true)
    return MVAResult(λ[order], V[:, order], B̄, size(B, 1), field)
end

function _moments(B, v::Val{N}) where {N}
    n = size(B, 1)
    B̄ = SVector{N}(sum(Bc) / n for Bc in eachcol(B))
    # Center before the products: the raw-moment form ⟨BBᵀ⟩ - B̄B̄ᵀ cancels catastrophically when |B̄| ≫ δB.
    M = zero(SMatrix{N, N, eltype(B̄)})
    @inbounds for k in axes(B, 1)
        d = SVector(ntuple(i -> B[k, i], v)) - B̄
        M += d * d'
    end
    return B̄, Symmetric(M / n)
end

function _moments(B, ::Nothing)
    B̄ = vec(sum(B; dims = 1)) / size(B, 1)
    C = B .- B̄'
    return B̄, Symmetric(C' * C / size(B, 1))
end

"""
    mva_eigen(B; dim = nothing, field = :B, check = false) -> F::MVAResult

Perform minimum variance analysis of the magnetic field `B` or maximum variance analysis of the electric field `E` when `field=:E`.

Samples lie along `dim`, by default the time dimension of a container (`DimArray`, `DataVariable`, …); plain matrices must pass `dim`.

`check=true` warns unless `issuccess(F)`.

## Notes

For a one-dimensional current layer, the tangential electric field components are approximately constant across the boundary, while the normal component exhibits the largest variation. Therefore, the eigenvector corresponding to the **maximum eigenvalue** ``λ_1`` (first column of `F.vectors`) gives an estimate of the boundary normal direction.
"""
function mva_eigen(B; dim = nothing, check = false, field = :B)
    field in (:B, :E) || throw(ArgumentError("field must be :B or :E, got $(repr(field))"))
    data = unwrap(B)
    F = _mva_eigen(resolve_dim(B, dim) == 1 ? data : transpose(data), field)
    check && !issuccess(F) && @warn "MVA$(F.field) eigenvalues are poorly separated; the normal is unreliable" F.values
    return F
end

# Unit extensions pass `B` stripped to `raw`, in units of `oneunit(first(B))`, so the numerics run on plain floats.
function _quantity_mva_eigen(B, raw, is_electric; field = nothing, kwargs...)
    scale = oneunit(first(B))
    field = @something field (is_electric ? :E : :B)
    F = mva_eigen(raw; field, kwargs...)
    return MVAResult(F.values .* scale^2, F.vectors, F.mean .* scale, F.n, F.field)
end

function issuccess(F::MVAResult; r0 = 5)
    _check3(F)
    λ = F.values
    return abs(F.field == :E ? λ[1] / λ[2] : λ[2] / λ[3]) > r0
end


"""
    Δφij(λᵢ, λⱼ, λ₃, M)

Calculate the phase error between components i and j according to:
|Δφᵢⱼ| = |Δφⱼᵢ| = √(λ₃/(M-1) * (λᵢ + λⱼ - λ₃)/(λᵢ - λⱼ)²)

With `λᵢ` eigenvalue and `M` number of samples.
"""
function Δφij(λᵢ, λⱼ, λ₃, M)
    return sqrt((λ₃ / (M - 1)) * (λᵢ + λⱼ - λ₃) / (λᵢ - λⱼ)^2)
end

"""
    normal_component_error(F::MVAResult)

Statistical error of the mean field's normal component, ``⟨𝐁⟩·𝐱₃`` for MVAB or ``⟨𝐄⟩·𝐱₁`` for MVAE:

``|Δ⟨𝐅·𝐱_n⟩| = \\sqrt{\\frac{λ_n}{M-1} + ∑_{k ≠ n} (Δφ_{nk} ⟨𝐅⟩·𝐱_k)^2}``
"""
function normal_component_error(F::MVAResult)
    _check3(F)
    λ, x, M = F.values, eachcol(F.vectors), F.n
    n = F.field == :E ? 1 : 3
    var = λ[n] / (M - 1)
    for k in 1:3
        k == n && continue
        var += (Δφij(λ[n], λ[k], λ[3], M) * dot(F.mean, x[k]))^2
    end
    return sqrt(var)
end

include("workload.jl")

end

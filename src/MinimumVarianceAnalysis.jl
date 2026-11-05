module MinimumVarianceAnalysis

using LinearAlgebra
import LinearAlgebra: issuccess
using StaticArrays
using TimeseriesUtilities: resolve_dim, unwrap, rebuild
export mva, mva_eigen, issuccess

struct MVAResult{V, M, B, N}
    values::V
    vectors::M
    mean::B
    n::Int
    field::Symbol
    noise::N
    noise_n::Int
end

MVAResult(values, vectors, mean, n, field) = MVAResult(values, vectors, mean, n, field, nothing, 0)

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
function _mva_eigen(B, field, background)
    return size(B, 2) == 3 ? _mva_eigen(B, Val(3), field, background) : _mva_eigen(B, nothing, field, background)
end

function _mva_eigen(B, N, field, background)
    B̄, M = _moments(B, N)
    C, nbg = _pooled_cov(background, N)
    λ, V = eigen(isnothing(C) ? M : Symmetric(M - C))
    order = sortperm(λ; rev = true)
    vectors = V[:, order]
    noise = isnothing(C) ? nothing : vectors' * C * vectors
    return MVAResult(λ[order], vectors, B̄, size(B, 1), field, noise, nbg)
end

_pooled_cov(::Nothing, N) = nothing, 0
function _pooled_cov(bgs, N)
    n = sum(b -> size(b, 1), bgs)
    return sum(b -> size(b, 1) * parent(_moments(b, N)[2]), bgs) / n, n
end

_samples(B, dim) = (data = unwrap(B); resolve_dim(B, dim) == 1 ? data : transpose(data))

_intervals(b::Union{Tuple, AbstractVector{<:AbstractArray}}) = b
_intervals(b) = (b,)

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
    mva_eigen(B; dim = nothing, field = :B, background = nothing, check = false) -> F::MVAResult

Perform minimum variance analysis of the magnetic field `B` or maximum variance analysis of the electric field `E` when `field=:E`.

`check=true` warns unless `issuccess(F)`.

## Notes

For a one-dimensional current layer, the tangential electric field components are approximately constant across the boundary, while the normal component exhibits the largest variation. Therefore, the eigenvector corresponding to the **maximum eigenvalue** ``λ_1`` (first column of `F.vectors`) gives an estimate of the boundary normal direction.
"""
function mva_eigen(B; dim = nothing, check = false, field = :B, background = nothing)
    field in (:B, :E) || throw(ArgumentError("field must be :B or :E, got $(repr(field))"))
    X = _samples(B, dim)
    bg = isnothing(background) ? nothing : map(b -> _samples(b, dim), _intervals(background))
    isnothing(bg) || all(b -> size(b, 2) == size(X, 2), bg) ||
        throw(DimensionMismatch("background must have the $(size(X, 2)) components of B"))
    F = _mva_eigen(X, field, bg)
    check && !issuccess(F) && @warn "MVA$(F.field) eigenvalues are poorly separated; the normal is unreliable" F.values
    return F
end

# Unit extensions pass `strip`, which converts a quantity to a plain number in units of `oneunit(first(B))`.
function _quantity_mva_eigen(B, strip, is_electric; field = nothing, background = nothing, kwargs...)
    scale = oneunit(first(B))
    field = @something field (is_electric ? :E : :B)
    bg = isnothing(background) ? nothing : map(b -> strip.(b), _intervals(background))
    F = mva_eigen(strip.(B); field, background = bg, kwargs...)
    noise = isnothing(F.noise) ? nothing : F.noise .* scale^2
    return MVAResult(F.values .* scale^2, F.vectors, F.mean .* scale, F.n, F.field, noise, F.noise_n)
end

# Signal variances along each eigenvector, noise covariance in the eigenbasis, and the inflation of noise×noise
# terms by background's own sampling error.
# Without a background, noise is assumed isotropic with variance λ₃.
function _signal_noise(F::MVAResult)
    λ = F.values
    isnothing(F.noise) && return λ .- λ[3], λ[3] * I, 1
    return max.(λ, zero(eltype(λ))), F.noise, 1 + F.n / F.noise_n
end

"""
    issuccess(F::MVAResult; r0 = 5)

Whether the normal is well determined: ``λ_i - λ_j > (r_0 - 1) V_j``, with ``(i, j) = (2, 3)`` for MVAB and ``(1, 2)``
for MVAE, and ``V_j`` the larger of the variances of `B` and of `background` along ``𝐱_j``. 
Without `background` this is ``λ_i/λ_j > r_0``.
"""
function issuccess(F::MVAResult; r0 = 5)
    _check3(F)
    λ = F.values
    s², Σ = _signal_noise(F)
    i, j = F.field == :E ? (1, 2) : (2, 3)
    return λ[i] - λ[j] > (r0 - 1) * (s²[j] + Σ[j, j])
end


"""
    normal_component_error(F::MVAResult; M = F.n)

Statistical error of the mean field's normal component, ``⟨𝐁⟩·𝐱₃`` for MVAB or ``⟨𝐄⟩·𝐱₁`` for MVAE:

``|Δ⟨𝐅·𝐱_n⟩| = \\sqrt{\\frac{s_n^2 + σ_n^2}{M-1} + ∑_{k ≠ n} (Δφ_{nk} ⟨𝐅⟩·𝐱_k)^2}``

with the angular uncertainty of eigenvector ``n`` toward ``k``

``|Δφ_{nk}| = \\sqrt{\\frac{σ_n^2 s_k^2 + σ_k^2 s_n^2 + (1 + n/n_b)(σ_n^2 σ_k^2 + σ_{nk}^4)}{(M-1)(λ_n - λ_k)^2}}``

with signal variance ``s_k^2`` along eigenvector ``k``, noise covariance ``σ_{nk}^2`` (``σ_k^2 ≡ σ_{kk}^2``) in the
eigenbasis, and ``1 + n/n_b`` the sampling error of a background of ``n_b`` samples against the ``n`` of `B`, assuming
equal correlation times. Without `background`, ``σ_{nk}^2 = λ_3 δ_{nk}``, ``n_b = ∞`` and
``s_k^2 = λ_k - λ_3``, which gives Sonnerup & Scheible (1998) ``|Δφ_{nk}| = \\sqrt{\\frac{λ_3}{M-1} \\frac{λ_n + λ_k - λ_3}{(λ_n - λ_k)^2}}``.
"""
function normal_component_error(F::MVAResult; M = F.n)
    _check3(F)
    λ, x = F.values, eachcol(F.vectors)
    s², Σ, f = _signal_noise(F)
    n = F.field == :E ? 1 : 3
    var = s²[n] + Σ[n, n]
    for k in 1:3
        k == n && continue
        Δφ² = (Σ[n, n] * s²[k] + Σ[k, k] * s²[n] + f * (Σ[n, n] * Σ[k, k] + Σ[n, k]^2)) / (λ[n] - λ[k])^2
        var += Δφ² * dot(F.mean, x[k])^2
    end
    return sqrt(var / (M - 1))
end

include("workload.jl")

end

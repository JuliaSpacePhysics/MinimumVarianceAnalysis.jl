# Maximum Variance Analysis of Electric Field (MVAE)
# Reference: Sonnerup, B. U. Ö., & Scheible, M. (1998). Minimum and maximum variance analysis.
#            ISSI Scientific Reports Series, 1, 185–220.

"""
    convection_efield(v, B; dim=nothing)

Compute the convection electric field ``\\mathbf{E} = -\\mathbf{v} × \\mathbf{B}`` from
plasma velocity `v` and magnetic field `B`.

This can be used as a proxy for the measured electric field when direct measurements are
unavailable (Sonnerup et al., 1987).

See also: [`mvae`](@ref), [`mvae_eigen`](@ref)
"""
function convection_efield(v, B; dim = nothing)
    dim = something(dim, 1)
    v_in = dim == 1 ? v : v'
    B_in = dim == 1 ? B : B'
    E = stack(-cross3(view(v_in, i, :), view(B_in, i, :)) for i in axes(v_in, 1); dims = 1)
    return dim == 1 ? E : E'
end

"""
    mvae_eigen(E::AbstractMatrix; dim=nothing, sort=(;), check=false) -> F::Eigen

Perform maximum variance analysis of the electric field `E` (Sonnerup & Scheible, 1998).

For a one-dimensional current layer, the tangential electric field components are
approximately constant across the boundary, while the normal component exhibits the
largest variation. Therefore, the eigenvector corresponding to the **maximum eigenvalue**
``λ_1`` (first column of `F.vectors`) gives an estimate of the boundary normal direction.

Return `Eigen` factorization object `F` which contains the eigenvalues in `F.values`
and the eigenvectors in the columns of the matrix `F.vectors`.

Set `check=true` to check the reliability of the result.

See also: [`mvae`](@ref), [`check_mvae_eigen`](@ref), [`convection_efield`](@ref)
"""
function mvae_eigen(E; dim = nothing, sort = (;), check = false)
    dim = something(dim, 1)
    in = dim == 1 ? E : E'
    N = size(in, 2)
    F = _mva_eigen(in, Val(N); sort)
    check && check_mvae_eigen(F)
    return F
end

function mvae_eigen(E::AbstractMatrix{Q}; kwargs...) where {Q <: Quantity}
    F = mvae_eigen(ustrip(E); kwargs...)
    return Eigen(F.values * unit(Q)^2, F.vectors)
end

"""
    mvae(V, E=V; dim=nothing, kwargs...)

Transform a timeseries `V` into the coordinate system obtained from maximum variance
analysis of electric field `E` along the `dim` dimension (time).

The columns of the result correspond to: maximum variance (normal), intermediate,
and minimum variance directions.

See also: [`mvae_eigen`](@ref), [`convection_efield`](@ref), [`transform`](@ref)
"""
function mvae(V, E = V; dim = nothing, kwargs...)
    F = mvae_eigen(E; dim, kwargs...)
    return transform(V, F; dim)
end

"""
    check_mvae_eigen(F; r0=5, verbose=false)

Check the quality of the MVAE result.

For MVAE, a reliable normal direction requires the maximum eigenvalue ``λ_1`` to be
well-separated from the intermediate eigenvalue ``λ_2``. The ratio ``|λ_1 / λ_2| > r_0``
is used as a quality indicator (default ``r_0 = 5``).

See also: [`mvae_eigen`](@ref)
"""
function check_mvae_eigen(F; r0 = 5, verbose = false)
    r = abs(F.values[1] / F.values[2])
    flag = r > r0
    verbose && begin
        println(F.vectors)
        println("Ratio of maximum variance to intermediate variance = ", r)
        flag && @info "Seems to be a proper MVAE attempt!"
    end
    flag || @warn "Take the MVAE result with a grain of salt!"
    return flag
end

################
# Error Estimate
################

"""
    E_x1_error(λ₁, λ₂, λ₃, M, E, x₁, x₂, x₃)

Calculate the composite statistical error estimate for ⟨E·x₁⟩ (the mean electric field
along the maximum variance / normal direction):

``|Δ⟨\\mathbf{E}·\\mathbf{x}_1⟩| = \\sqrt{\\frac{λ_1}{M-1} + (Δφ_{12}⟨\\mathbf{E}⟩·\\mathbf{x}_2)^2 + (Δφ_{13}⟨\\mathbf{E}⟩·\\mathbf{x}_3)^2}``

Parameters:

  - λ₁, λ₂, λ₃: eigenvalues in descending order
  - M: number of samples
  - E: mean electric field vector
  - x₁, x₂, x₃: eigenvectors
"""
function E_x1_error(λ₁, λ₂, λ₃, M, E, x₁, x₂, x₃)
    Δφ₁₂ = Δφij(λ₁, λ₂, λ₃, M)
    Δφ₁₃ = Δφij(λ₁, λ₃, λ₃, M)
    E_x₂ = dot(E, x₂)
    E_x₃ = dot(E, x₃)
    return sqrt(λ₁ / (M - 1) + (Δφ₁₂ * E_x₂)^2 + (Δφ₁₃ * E_x₃)^2)
end

E_x1_error(F::Eigen, M, E) = E_x1_error(F.values..., M, E, eachcol(F.vectors)...)

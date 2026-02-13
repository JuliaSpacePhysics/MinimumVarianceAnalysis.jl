# Maximum Variance Analysis of Electric Field (MVAE)

"""
    convection_efield(v, B; dim=nothing)

Compute the convection electric field ``\\mathbf{E} = -\\mathbf{v} × \\mathbf{B}`` from
plasma velocity `v` and magnetic field `B`.

This can be used as a proxy for the measured electric field when direct measurements are
unavailable.
"""
function convection_efield(v, B; dim = nothing)
    dim = something(dim, 1)
    v_in = dim == 1 ? v : v'
    B_in = dim == 1 ? B : B'
    E = stack(-cross3(view(v_in, i, :), view(B_in, i, :)) for i in axes(v_in, 1); dims = 1)
    return dim == 1 ? E : E'
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

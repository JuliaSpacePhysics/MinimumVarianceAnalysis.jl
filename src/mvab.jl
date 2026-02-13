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

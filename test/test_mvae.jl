@testitem "convection_efield" begin
    using LinearAlgebra

    # v = [1,0,0], B = [0,1,0] → E = -v×B = -(0,0,1) = [0,0,-1]
    # v = [0,1,0], B = [1,0,0] → E = -v×B = -(0,0,-1) = [0,0,1]
    v = [1.0 0.0 0.0; 0.0 1.0 0.0]
    B = [0.0 1.0 0.0; 1.0 0.0 0.0]
    E = convection_efield(v, B)
    @test E ≈ [0.0 0.0 -1.0; 0.0 0.0 1.0]

    # Test with dim=2
    E2 = convection_efield(collect(v'), collect(B'); dim = 2)
    @test E2 ≈ collect(E')
end

@testitem "E_x1_error" begin
    using LinearAlgebra
    using MinimumVarianceAnalysis: E_x1_error, Δφij

    λ₁, λ₂, λ₃ = 10.0, 2.0, 0.5
    M = 100
    E = [1.0, 0.5, 0.1]
    x₁ = [1.0, 0.0, 0.0]
    x₂ = [0.0, 1.0, 0.0]
    x₃ = [0.0, 0.0, 1.0]

    err = E_x1_error(λ₁, λ₂, λ₃, M, E, x₁, x₂, x₃)
    @test err > 0
    @test isfinite(err)

    # Test convenience method with Eigen
    F = Eigen([λ₁, λ₂, λ₃], hcat(x₁, x₂, x₃))
    err2 = E_x1_error(F, M, E)
    @test err ≈ err2

    # Error should increase with fewer samples
    err_fewer = E_x1_error(λ₁, λ₂, λ₃, 10, E, x₁, x₂, x₃)
    @test err_fewer > err
end
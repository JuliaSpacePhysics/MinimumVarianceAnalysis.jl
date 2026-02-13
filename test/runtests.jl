using TestItemRunner

@run_package_tests

@testitem "Aqua" begin
    using Aqua
    Aqua.test_all(MinimumVarianceAnalysis)
end

@testitem "mva_eigen basic functionality with synthetic data" begin
    using LinearAlgebra
    using DimensionalData
    using Unitful

    # Test with a simple 3D magnetic field data where we know the expected result
    # Create synthetic data with known variance structure
    Bl = cosd.(0:30:180)
    Bm = sind.(0:30:180)
    Bn = 0.1 .* ones(7)
    B = hcat(Bl, Bm, Bn)

    F = mva_eigen(B; check = true)

    # Check that we get 3 eigenvalues and 3 eigenvectors
    @test length(F.values) == 3
    @test size(F.vectors) == (3, 3)
    # Check that eigenvalues are sorted in descending order (default behavior)
    @test F.values[1] >= F.values[2] >= F.values[3]
    @test F.vectors[:, 1] ≈ [1.0, 0.0, 0.0] || F.vectors[:, 1] ≈ [-1.0, 0.0, 0.0]
    @test F.vectors[:, 2] ≈ [0.0, 1.0, 0.0] || F.vectors[:, 2] ≈ [0.0, -1.0, 0.0]
    @test F.vectors[:, 3] ≈ [0.0, 0.0, 1.0] || F.vectors[:, 3] ≈ [0.0, 0.0, -1.0]
    # Check that eigenvectors are orthonormal
    @test F.vectors' * F.vectors ≈ I
    @test abs(det(F.vectors)) ≈ 1
    @test mva(B) .* diag(F.vectors)' ≈ B # Note: eigenvectors are only unique up to sign

    # Dimensional data
    B1 = DimArray(B * u"nT", (Ti(0:6), Y(1:3)))
    B2 = DimArray(collect(B') * u"nT", (Y(1:3), Ti(0:6)))
    F1 = mva_eigen(B1)
    F2 = mva_eigen(B2)
    @test F1.values ≈ F2.values
    @test F1.vectors ≈ F2.vectors

    # Test rotation
    B1_rot = mva(B1)
    B2_rot = mva(B2)
    @test B1_rot == B2_rot'

    # Test edge cases and robustness
    B_const = ones(10, 3)
    F_const = mva_eigen(B_const)
    @test F_const.values == [0, 0, 0]
end

@testitem "mvae_eigen basic functionality with synthetic data" begin
    using LinearAlgebra
    using DimensionalData
    using Unitful

    # Create synthetic E-field data with known variance structure
    # Normal component (max variance) along z, tangential components (low variance) along x, y
    En = cosd.(0:30:180)
    Et1 = 0.1 .* ones(7)
    Et2 = 0.05 .* ones(7)
    E = hcat(Et1, Et2, En)

    F = mvae_eigen(E; check = true)

    # Check that we get 3 eigenvalues and 3 eigenvectors
    @test length(F.values) == 3
    @test size(F.vectors) == (3, 3)
    # Check that eigenvalues are sorted in descending order (default behavior)
    @test F.values[1] >= F.values[2] >= F.values[3]
    # Maximum variance direction (column 1) should be the z-axis (normal)
    @test F.vectors[:, 1] ≈ [0.0, 0.0, 1.0] || F.vectors[:, 1] ≈ [0.0, 0.0, -1.0]
    # Check that eigenvectors are orthonormal
    @test F.vectors' * F.vectors ≈ I
    @test abs(det(F.vectors)) ≈ 1
    # Check transformation
    @test mvae(E) .* diag(F.vectors)' ≈ E

    # Dimensional data
    E1 = DimArray(E * u"mV/m", (Ti(0:6), Y(1:3)))
    E2 = DimArray(collect(E') * u"mV/m", (Y(1:3), Ti(0:6)))
    F1 = mvae_eigen(E1)
    F2 = mvae_eigen(E2)
    @test F1.values ≈ F2.values
    @test F1.vectors ≈ F2.vectors

    # Test with dim=2
    E_t = collect(E')
    Ft = mvae_eigen(E_t; dim = 2)
    @test Ft.values ≈ F.values
    @test Ft.vectors ≈ F.vectors

    # Test edge case: constant field
    E_const = ones(10, 3)
    F_const = mvae_eigen(E_const)
    @test F_const.values == [0, 0, 0]
end

@testitem "check_mvae_eigen quality check" begin
    using LinearAlgebra: Eigen

    # Good MVAE: well-separated maximum eigenvalue
    F_good = Eigen([100.0, 5.0, 1.0], [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
    @test check_mvae_eigen(F_good; r0 = 5) == true

    # Bad MVAE: poorly separated eigenvalues
    F_bad = Eigen([3.0, 2.5, 1.0], [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
    @test check_mvae_eigen(F_bad; r0 = 5) == false
end

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

    # Test that E = -v × B satisfies the identity: v ⋅ E = 0
    v3 = randn(20, 3)
    B3 = randn(20, 3)
    E3 = convection_efield(v3, B3)
    for i in axes(v3, 1)
        @test dot(view(v3, i, :), view(E3, i, :)) ≈ 0.0 atol = 1e-12
    end
end

@testitem "mvae with convection_efield integration" begin
    using LinearAlgebra

    # Create a scenario where the normal direction (z) has large E variation
    # and tangential directions have small E variation
    n = 50
    t = range(0, π, length = n)

    # Plasma velocity: mostly in x-direction, varying magnitude
    v = hcat(10.0 .* cos.(t), 0.5 .* ones(n), 0.2 .* ones(n))

    # Magnetic field: mostly in y-direction
    B = hcat(0.1 .* ones(n), 5.0 .* sin.(t), 0.3 .* ones(n))

    # Compute convection E-field
    E = convection_efield(v, B)
    @test size(E) == (n, 3)

    # Should be able to run MVAE on this
    F = mvae_eigen(E)
    @test length(F.values) == 3
    @test F.values[1] >= F.values[2] >= F.values[3]
    @test F.vectors' * F.vectors ≈ I

    # Transform should preserve dimensions
    E_transformed = mvae(E)
    @test size(E_transformed) == size(E)
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

@testitem "right-handed coordinate system check" begin
    using LinearAlgebra: Eigen
    using MinimumVarianceAnalysis: is_right_handed

    v1 = [1.0, 0.0, 0.0]
    v2 = [0.0, 1.0, 0.0]
    v3 = [0.0, 0.0, 1.0]
    v3_left = [0.0, 0.0, -1.0]

    # Test with Eigen object with right-handed system
    F_right = Eigen([3.0, 2.0, 1.0], [v1 v2 v3])
    @test is_right_handed(F_right) == true

    # Test with Eigen object with left-handed system
    F_left = Eigen([3.0, 2.0, 1.0], [v1 v2 v3_left])
    @test is_right_handed(F_left) == false
end
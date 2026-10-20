using TestItemRunner

@run_package_tests

@testitem "Aqua" begin
    using Aqua
    Aqua.test_all(MinimumVarianceAnalysis)
end

@testitem "MVAB basic functionality with synthetic data" begin
    using LinearAlgebra
    using DimensionalData
    using Unitful

    # Test with a simple 3D magnetic field data where we know the expected result
    # Create synthetic data with known variance structure
    Bl = cosd.(0:30:180)
    Bm = sind.(0:30:180)
    Bn = 0.1 .* ones(7)
    B = hcat(Bl, Bm, Bn)

    F = mva_eigen(B; dim = 1, check = true)

    @test F.vectors[:, 1] ≈ [1.0, 0.0, 0.0] || F.vectors[:, 1] ≈ [-1.0, 0.0, 0.0]
    @test F.vectors[:, 2] ≈ [0.0, 1.0, 0.0] || F.vectors[:, 2] ≈ [0.0, -1.0, 0.0]
    @test F.vectors[:, 3] ≈ [0.0, 0.0, 1.0] || F.vectors[:, 3] ≈ [0.0, 0.0, -1.0]
    @test F.mean ≈ vec(sum(B; dims = 1)) / 7
    @test mva(B; dim = 1) .* diag(F.vectors)' ≈ B # Note: eigenvectors are only unique up to sign

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
    @test B1_rot ≈ B2_rot'
    @test dims(B1_rot, Ti) == dims(B1, Ti)
    @test !hasdim(B1_rot, Y) # components lie along the eigenvectors now
    @test mva(B2, F1) ≈ B2_rot

    @test_throws ArgumentError mva_eigen(B)

    # Test edge cases and robustness
    B_const = ones(10, 3)
    F_const = mva_eigen(B_const; dim = 1)
end

@testitem "MVAE basic functionality with synthetic data" begin
    using LinearAlgebra
    using DimensionalData
    using Unitful

    # Create synthetic E-field data with known variance structure
    # Normal component (max variance) along z, tangential components (low variance) along x, y
    En = cosd.(0:30:180)
    Et1 = 0.1 .* ones(7)
    Et2 = 0.05 .* ones(7)
    E = hcat(Et1, Et2, En)

    F = mva_eigen(E; dim = 1, check = true, field = :E)

    # Check that eigenvalues are sorted in descending order (default behavior)
    @test F.values[1] >= F.values[2] >= F.values[3]
    # Maximum variance direction (column 1) should be the z-axis (normal)
    @test F.vectors[:, 1] ≈ [0.0, 0.0, 1.0] || F.vectors[:, 1] ≈ [0.0, 0.0, -1.0]

    # Dimensional data
    E1 = DimArray(E * u"mV/m", (Ti(0:6), Y(1:3)))
    E2 = DimArray(collect(E') * u"mV/m", (Y(1:3), Ti(0:6)))
    F1 = mva_eigen(E1; field = :E)
    F2 = mva_eigen(E2; field = :E)
    @test F1.values ≈ F2.values
    @test F1.vectors ≈ F2.vectors

    # Test with dim=2
    E_t = collect(E')
    Ft = mva_eigen(E_t; dim = 2, field = :E)
    @test Ft.values ≈ F.values
    @test Ft.vectors ≈ F.vectors
end

@testitem "issuccess quality ratio" begin
    using LinearAlgebra: I
    using MinimumVarianceAnalysis: MVAResult

    F(λ, field) = MVAResult(λ, Matrix(1.0I, length(λ), length(λ)), zeros(length(λ)), 10, field)
    # MVAE needs λ₁ ≫ λ₂; MVAB needs λ₂ ≫ λ₃
    @test issuccess(F([100.0, 5.0, 1.0], :E))
    @test !issuccess(F([100.0, 5.0, 1.0], :B))
    @test issuccess(F([3.0, 2.5, 0.1], :B))
    @test !issuccess(F([3.0, 2.5, 0.1], :E))
    @test issuccess(F([3.0, 2.5, 0.1], :E); r0 = 1.1)
    # λ₃ must be the smallest eigenvalue, not the third of four
    @test_throws ArgumentError issuccess(F([4.0, 3.0, 2.0, 1.0], :B))
end

@testitem "normal and normal_component_error follow field" begin
    using LinearAlgebra: I
    using MinimumVarianceAnalysis: MVAResult, normal, normal_component_error, Δφij

    λ₁, λ₂, λ₃, M = 10.0, 2.0, 0.5, 100
    m = [1.0, 0.5, 0.1]
    FB = MVAResult([λ₁, λ₂, λ₃], Matrix(1.0I, 3, 3), m, M, :B)
    FE = MVAResult([λ₁, λ₂, λ₃], Matrix(1.0I, 3, 3), m, M, :E)
    @test normal(FB) == [0, 0, 1]
    @test normal(FE) == [1, 0, 0]
    # error of ⟨B·x₃⟩ (MVAB) and ⟨E·x₁⟩ (MVAE)
    @test normal_component_error(FB) ≈ sqrt(λ₃ / (M - 1) + (Δφij(λ₃, λ₂, λ₃, M) * m[2])^2 + (Δφij(λ₃, λ₁, λ₃, M) * m[1])^2)
    @test normal_component_error(FE) ≈ sqrt(λ₁ / (M - 1) + (Δφij(λ₁, λ₂, λ₃, M) * m[2])^2 + (Δφij(λ₁, λ₃, λ₃, M) * m[3])^2)
end

@testitem "right-handed coordinate system check" begin
    using MinimumVarianceAnalysis: MVAResult, is_right_handed

    F(v3) = MVAResult([3.0, 2.0, 1.0], [[1.0, 0.0, 0.0] [0.0, 1.0, 0.0] v3], zeros(3), 10, :B)
    @test is_right_handed(F([0.0, 0.0, 1.0]))
    @test !is_right_handed(F([0.0, 0.0, -1.0]))
end

@testitem "more than 3 components" begin
    using LinearAlgebra
    using MinimumVarianceAnalysis: _moments

    B = [cosd.(0:30:180) sind.(0:30:180) 0.1 .* ones(7)]
    @test all(map(≈, _moments(B, nothing), _moments(B, Val(3))))

    B4 = [3 .* cosd.(0:30:330) 2 .* sind.(0:30:330) 0.5 .* cosd.(0:60:660) 0.1 .* ones(12)]
    F = mva_eigen(B4; dim = 1)
    @test abs.(F.vectors) ≈ I
end

@testitem "containers resolve their time dimension" begin
    using TimeseriesUtilities: Timed, unwrap

    B = [cosd.(0:30:180) sind.(0:30:180) 0.1 .* ones(7)]
    F = mva_eigen(B; dim = 1)
    Ft = mva_eigen(Timed(collect(B'), 0:6; dim = 2))
    @test Ft.values ≈ F.values
    @test Ft.vectors ≈ F.vectors
    @test unwrap(mva(Timed(collect(B'), 0:6; dim = 2))) ≈ mva(B; dim = 1)'
end

@testitem "large mean field in Float32" begin
    using LinearAlgebra
    using MinimumVarianceAnalysis: _moments

    # |B̄| ≫ δB as at a magnetopause crossing; raw second moments lose λ₃ to Float32 cancellation
    B = [40 .+ cosd.(0:719) -30 .+ 0.3 .* sind.(0:719) 20 .+ 0.03 .* cosd.(0:3:2157)]
    λ₃ = mva_eigen(B; dim = 1).values[3]
    @test mva_eigen(Float32.(B); dim = 1).values[3] ≈ λ₃ rtol = 1.0e-4
    @test eigmin(_moments(Float32.(B), nothing)[2]) ≈ λ₃ rtol = 1.0e-4
end

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
    using MinimumVarianceAnalysis: MVAResult, normal, normal_component_error

    λ₁, λ₂, λ₃, M = 10.0, 2.0, 0.5, 100
    m = [1.0, 0.5, 0.1]
    FB = MVAResult([λ₁, λ₂, λ₃], Matrix(1.0I, 3, 3), m, M, :B)
    FE = MVAResult([λ₁, λ₂, λ₃], Matrix(1.0I, 3, 3), m, M, :E)
    @test normal(FB) == [0, 0, 1]
    @test normal(FE) == [1, 0, 0]
    # Sonnerup & Scheible (1998) Eqs. 8.23-8.24 for ⟨B·x₃⟩ (MVAB) and ⟨E·x₁⟩ (MVAE)
    Δφ(λᵢ, λⱼ) = sqrt(λ₃ / (M - 1) * (λᵢ + λⱼ - λ₃) / (λᵢ - λⱼ)^2)
    @test normal_component_error(FB) ≈ sqrt(λ₃ / (M - 1) + (Δφ(λ₃, λ₂) * m[2])^2 + (Δφ(λ₃, λ₁) * m[1])^2)
    @test normal_component_error(FE) ≈ sqrt(λ₁ / (M - 1) + (Δφ(λ₁, λ₂) * m[2])^2 + (Δφ(λ₁, λ₃) * m[3])^2)
    @test normal_component_error(FB; M = 10) ≈ normal_component_error(FB) * sqrt((M - 1) / 9)
end

@testitem "background removes the bias of anisotropic noise" begin
    using LinearAlgebra
    using Random: Xoshiro
    using Unitful
    using MinimumVarianceAnalysis: normal, normal_component_error

    # Normal along z; noise, like Alfvénic turbulence, has almost no variance along e, 45° from the normal
    rng = Xoshiro(1)
    e = normalize([0.0, 1.0, 1.0])
    Σ½ = 0.3 * (I - e * e') + 0.01 * e * e'
    noise(N) = randn(rng, N, 3) * Σ½
    t = range(-3, 3; length = 20_000)
    B = [tanh.(t) 0.3 .* sech.(t) fill(0.5, length(t))] + noise(length(t))
    before, after = noise(20_000) .+ [-1 0 0.5], noise(20_000) .+ [1 0 0.5]
    angle(n) = acosd(abs(n[3]))

    F0 = mva_eigen(B; dim = 1)
    # wrong, yet reported reliable
    @test angle(normal(F0)) > 20 && issuccess(F0)
    @test abs(abs(F0.mean ⋅ normal(F0)) - 0.5) > 3normal_component_error(F0)
    F = mva_eigen(B; dim = 1, background = (before, after))
    @test angle(normal(F)) < 3
    @test F.noise[3, 3] ≈ normal(F)' * Σ½^2 * normal(F) rtol = 0.05

    # a background in other units is converted to those of B
    Fu = mva_eigen(B .* u"nT"; dim = 1, background = (before .* 1.0e-9u"T", after .* u"nT"))
    @test ustrip.(u"nT^2", Fu.values) ≈ F.values
    @test ustrip(u"nT", normal_component_error(Fu)) ≈ normal_component_error(F)
end

@testitem "normal_component_error matches the scatter under anisotropic noise and a short background" begin
    using LinearAlgebra
    using Random: Xoshiro
    using Statistics: std, median
    using MinimumVarianceAnalysis: normal, normal_component_error

    # Both the noise cross-covariance between eigenvectors and the background's sampling error each add ~30% here.
    rng = Xoshiro(2)
    e = normalize([0.0, 1.0, 1.0])
    Σ½ = 0.3 * (I - e * e') + 0.01 * e * e'
    noise(N) = randn(rng, N, 3) * Σ½
    t = range(-3, 3; length = 2000)
    S = [tanh.(t) 0.3 .* sech.(t) fill(0.5, length(t))]
    Fs = [mva_eigen(S + noise(length(t)); dim = 1, background = (noise(1000) .+ [-1 0 0.5], noise(1000) .+ [1 0 0.5])) for _ in 1:200]
    Bn = [abs(F.mean ⋅ normal(F)) for F in Fs]
    @test std(Bn) / median(normal_component_error.(Fs)) ≈ 1 atol = 0.2
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

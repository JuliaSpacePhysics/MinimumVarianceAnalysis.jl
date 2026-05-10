@testitem "Unitful input infers field" begin
    using Unitful

    B = hcat(cosd.(0:30:180), sind.(0:30:180), 0.1 .* ones(7))
    F = mva_eigen(B .* u"T")

    @test Unitful.ustrip.(F.values) ≈ mva_eigen(B).values
    @test all(Unitful.unit.(F.values) .== Ref(u"T^2"))
    @test F.vectors ≈ mva_eigen(B).vectors

    E = hcat(0.1 .* ones(7), 0.05 .* ones(7), cosd.(0:30:180))
    FE = mva_eigen(E .* u"mV/m")
    @test FE.vectors[:, 1] ≈ [0.0, 0.0, 1.0] || FE.vectors[:, 1] ≈ [0.0, 0.0, -1.0]
end

@testitem "Quantity packages infer field and preserve units" begin
    using DynamicQuantities
    using FlexUnits
    import FlexUnits.UnitRegistry

    B = hcat(cosd.(0:30:180), sind.(0:30:180), 0.1 .* ones(7))
    F = mva_eigen(B)

    dq_B = B .* DynamicQuantities.uparse("T")
    dq_F = mva_eigen(dq_B)
    @test DynamicQuantities.ustrip.(dq_F.values) ≈ F.values
    @test all(DynamicQuantities.dimension.(dq_F.values) .== Ref(DynamicQuantities.dimension(DynamicQuantities.uparse("T")^2)))
    @test dq_F.vectors ≈ F.vectors

    flex_B = B .* UnitRegistry.uparse("T")
    flex_F = mva_eigen(flex_B)
    @test FlexUnits.ustrip.(flex_F.values) ≈ F.values
    @test all(FlexUnits.dimension.(flex_F.values) .== Ref(FlexUnits.dimension(UnitRegistry.uparse("T")^2)))
    @test flex_F.vectors ≈ F.vectors

    E = hcat(0.1 .* ones(7), 0.05 .* ones(7), cosd.(0:30:180))
    dq_E = mva_eigen(E .* DynamicQuantities.uparse("V/m"))
    flex_E = mva_eigen(E .* UnitRegistry.uparse("V/m"))
    @test dq_E.vectors[:, 1] ≈ [0.0, 0.0, 1.0] || dq_E.vectors[:, 1] ≈ [0.0, 0.0, -1.0]
    @test flex_E.vectors[:, 1] ≈ [0.0, 0.0, 1.0] || flex_E.vectors[:, 1] ≈ [0.0, 0.0, -1.0]
end

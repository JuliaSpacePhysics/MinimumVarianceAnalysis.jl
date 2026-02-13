@testitem "Unitful input infers field" begin
    using MinimumVarianceAnalysis: _field
    using Unitful
    @test _field(u"mV/m") == :E
    @test _field(u"T") == :B
end
using MinimumVarianceAnalysis
using Documenter

makedocs(;
    modules = [MinimumVarianceAnalysis],
    sitename = "MinimumVarianceAnalysis.jl",
    format = Documenter.HTML(;),
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo = "github.com/JuliaSpacePhysics/MinimumVarianceAnalysis.jl"
)

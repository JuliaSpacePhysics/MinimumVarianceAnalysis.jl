using MinimumVarianceAnalysis
using Documenter

DocMeta.setdocmeta!(MinimumVarianceAnalysis, :DocTestSetup, :(using MinimumVarianceAnalysis); recursive=true)

makedocs(;
    modules=[MinimumVarianceAnalysis],
    authors="Beforerr <zzj956959688@gmail.com> and contributors",
    sitename="MinimumVarianceAnalysis.jl",
    format=Documenter.HTML(;
        canonical="https://Beforerr.github.io/MinimumVarianceAnalysis.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Beforerr/MinimumVarianceAnalysis.jl",
    devbranch="main",
)

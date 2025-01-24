using Documenter, Fliwer

makedocs(
    sitename = "Fliwer.jl",
    pages = [
        "index.md",
        "Gradient Test" => "gradient.md",
        "Poisson 2D - 1 Phase" => "poisson.md",
        "Poisson 2D - 2 Phases" => "poisson_2phases.md",
        "Heat 2D - 1 Phase - Robin" => "heat_2D_Robin.md",
        "Heat 1D - 2 Phases" => "heat_1D_2ph.md",
    ]
)

deploydocs(
    repo="github.com/Fastaxx/Fliwer.jl.git",
    versions = nothing
)
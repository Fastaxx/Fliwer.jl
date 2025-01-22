using Documenter, Fliwer

makedocs(
    sitename = "Fliwer.jl",
    pages = [
        "index.md",
        "Gradient Test" => "gradient.md",
        "Poisson 2D - 1 Phase" => "poisson.md"
    ]
)

deploydocs(
    repo="github.com/Fastaxx/Fliwer.jl.git"
)
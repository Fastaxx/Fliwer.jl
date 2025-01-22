using Documenter, Fliwer

makedocs(
    sitename = "Fliwer.jl",
    pages = [
        "index.md",
        "Poisson 2D - 1 Phase" => "poisson.md"
    ]
)

deploydocs(
    repo="github.com/Fastaxx/Fliwer.jl.git"
)
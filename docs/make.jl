using Documenter, Fliwer

makedocs(sitename="Fliwer Documentation", authors = "Louis Libat", pages= ["index.md", "Examples" => "examples.md"])

deploydocs(repo = "github.com/Fastaxx/Fliwer.jl.git", target = "build", branch = "gh-pages", push_preview = false)
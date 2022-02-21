using Documenter
using Revise
using Smg2s
using SparseArrays

makedocs(
        modules = [Smg2s],
        sitename="Smg2s.jl",
        format = Documenter.HTML(
              prettyurls = get(ENV, "CI", nothing) == "true",
        ),
        pages = [
                "Home" => "index.md",
                "User Guide" => Any[
                        "Getting Started" => "getting_started.md"
                        "Customization" => "custorm.md"
                        "Examples" => "examples.md"
                ],
                "Gallery: Sparsity Patterns" => "gallery.md",
                "Citing Smg2s" => "citing.md",
                "Licence" => "licence.md"
        ]


)

deploydocs(
    repo = "https://github.com/Smg2s/Smg2s.jl")

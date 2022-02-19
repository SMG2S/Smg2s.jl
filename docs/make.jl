using Documenter
using Revise
using SMG2S
using SparseArrays

makedocs(
        modules = [SMG2S],
        sitename="SMG2S.jl",
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
                "Citing SMG2S" => "citing.md",
                "Licence" => "licence.md"
        ]


)

deploydocs(
    repo = "https://github.com/SMG2S/SMG2S.jl")

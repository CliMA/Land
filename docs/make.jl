using Documenter
using StomataModels

pages = Any[
    "Home" => "index.md",
    "API"  => "API.md"
    ]


mathengine = MathJax(Dict(
    :TeX => Dict(
        :equationNumbers => Dict(:autoNumber => "AMS"),
        :Macros => Dict(),
    ),
))

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    mathengine = mathengine,
    collapselevel = 1,
)

makedocs(
    sitename = "StomataModels",
    format = format,

    clean = false,
    modules = [StomataModels],
    pages = pages,
)

deploydocs(
    repo = "github.com/Yujie-W/StomataModels.jl.git",
    target = "build",
    push_preview = true,
)
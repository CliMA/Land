using Documenter
using PlantHydraulics

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
    sitename = "PlantHydraulics",
    format = format,
    clean = false,
    modules = [PlantHydraulics],
    pages = pages,
)

deploydocs(
    repo = "github.com/Yujie-W/PlantHydraulics.jl.git",
    target = "build",
    devbranch = "main",
    push_preview = true,
)

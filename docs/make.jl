using Documenter
using WaterPhysics

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
    sitename = "WaterPhysics",
    format = format,
    clean = false,
    modules = [WaterPhysics],
    pages = pages,
)

deploydocs(
    repo = "github.com/Yujie-W/WaterPhysics.jl.git",
    target = "build",
    push_preview = true,
)

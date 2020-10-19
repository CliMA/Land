using Documenter
using CanopyLayers

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
    sitename = "CanopyLayers",
    format = format,
    clean = false,
    modules = [CanopyLayers],
    pages = pages,
)

deploydocs(
    repo = "github.com/Yujie-W/CanopyLayers.jl.git",
    target = "build",
    devbranch = "main",
    push_preview = true,
)

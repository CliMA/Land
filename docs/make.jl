using Documenter
using Plants

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
    sitename = "Plants",
    format = format,

    clean = false,
    modules = [Plants],
    pages = pages,
)

deploydocs(
    repo = "github.com/Yujie-W/Plants.jl.git",
    target = "build",
    push_preview = true,
)
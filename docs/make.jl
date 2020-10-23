using Documenter
using Photosynthesis

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
    sitename = "Photosynthesis",
    format = format,
    clean = false,
    modules = [Photosynthesis],
    pages = pages,
)

deploydocs(
    repo = "github.com/Yujie-W/Photosynthesis.jl.git",
    target = "build",
    devbranch = "main",
    push_preview = true,
)

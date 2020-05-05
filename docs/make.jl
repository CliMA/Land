using Land, Documenter

pages = Any[
    "Home" => "index.md",
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
    sitename = "Land",
    format = format,
    clean = true,
    modules = [Documenter, Land],
    pages = pages,
)

deploydocs(
    repo = "github.com/climate-machine/Land.git",
    target = "build",
    push_preview = true,
)

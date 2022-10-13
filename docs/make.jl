using Documenter
using Land


# define default docs pages
pages = Any[
    "Home" => "index.md",
    "APIs" => "API.md"  ,
    "Tips" => "tips.md" ,
]


# format the docs
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
    assets = ["assets/favicon.ico"]
)


# build the docs
makedocs(
    sitename = "Land",
    format = format,
    clean = false,
    modules = [Land],
    pages = pages,
)


# deploy the docs to Github gh-pages
deploydocs(
    repo = "github.com/CliMA/Land.git",
    target = "build",
    devbranch = "main",
    push_preview = true,
)

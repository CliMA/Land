using Land, Documenter, Literate

generated_dir = joinpath(@__DIR__, "src", "generated") # generated files directory
rm(generated_dir, force = true, recursive = true)
mkpath(generated_dir)

include("list_of_tutorials.jl")          # defines a dict `tutorials`

pages = Any[
    "Home"           => "index.md",
    "Tips"           => "pages/Tips.md",
    "Hydraulics"     => "pages/Hydraulics.md",
    "Photosynthesis" => "pages/Photosynthesis.md",
    "CanopyRT"       => "pages/CanopyRT.md",
    "Utils"          => "pages/Utils.md",
    "Tutorials"      => tutorials,
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

    clean = false,
    #modules = [Documenter, Land],
    modules = [Land],    # remove Documenter to avoid tons of WARNINGS associated
    pages = pages,
)

deploydocs(
    repo = "github.com/CliMA/Land.git",
    target = "build",
    push_preview = true,
)

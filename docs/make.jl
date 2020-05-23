using Land, Documenter, Literate

generated_dir = joinpath(@__DIR__, "src", "generated") # generated files directory
rm(generated_dir, force = true, recursive = true)
mkpath(generated_dir)

include("list_of_tutorials.jl")          # defines a dict `tutorials`

pages = Any[
    "Home" => "index.md",
    "Plant Module" => "plant_module.md",
    "PhotosynthesisModels Module" => "photosynthesis_models.md",
    "Tutorials" => tutorials,
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
    modules = [Documenter, Land],
    pages = pages,
)

deploydocs(
    repo = "github.com/CliMA/Land.git",
    target = "build",
    push_preview = true,
)

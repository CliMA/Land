using Documenter
using Land




# define default docs pages
pages = Any[
    "Home" => "index.md",
    "Tips" => "tips.md" ,
    "APIs" => [
              "CanopyLayers"          => "CanopyLayers.md"         ,
              "Photosynthesis"        => "Photosynthesis.md"       ,
              "PlantHydraulics"       => "PlantHydraulics.md"      ,
              "StomataModels"         => "StomataModels.md"        ,
              "SoilPlantAirContinuum" => "SoilPlantAirContinuum.md",
              "Land"                  => "Land.md"
    ]
]




# add example and tutorial pages
# currently unavailable
generated_dir = joinpath(@__DIR__, "src", "generated") # generated files directory
rm(generated_dir, force = true, recursive = true)
mkpath(generated_dir)
#include("list_of_tutorials.jl")          # defines a dict `tutorials`




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

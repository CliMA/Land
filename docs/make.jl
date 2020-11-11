using Documenter
using Land




# define default docs pages
pages = Any[
    "Home" => "index.md",
    "Tips" => "tips.md" ,
    "APIs" => [
              "CanopyLayers"    => "submodules/CanopyLayers.md"         ,
              "Photosynthesis"  => "submodules/Photosynthesis.md"       ,
              "PlantHydraulics" => "submodules/PlantHydraulics.md"      ,
              "StomataModels"   => "submodules/StomataModels.md"        ,
              "SPAC"            => "submodules/SoilPlantAirContinuum.md",
              "Land"            => "submodules/Land.md"                 ,
    ],
]




# add example pages
generated_dir = joinpath(@__DIR__, "src", "generated")
rm(generated_dir, force = true, recursive = true)
mkpath(generated_dir)
include("examples.jl")




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

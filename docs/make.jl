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
];


# format the docs
mathengine = MathJax(Dict(
    :TeX => Dict(
        :equationNumbers => Dict(:autoNumber => "AMS"),
        :Macros => Dict(),
    ),
));

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    mathengine = mathengine,
    collapselevel = 1,
    assets = ["assets/favicon.ico"]
);


# build the docs
makedocs(
    sitename = "clima-land-v0.1",
    format = format,
    clean = false,
    # modules = [Land],
    pages = pages,
)


# deploy the docs to Github gh-pages
deploydocs(
    repo = "github.com/silicormosia/clima-land-v0.1.git",
    target = "build",
    devbranch = "main",
    push_preview = true,
)

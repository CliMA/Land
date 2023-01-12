using Documenter
using Land


# define default docs pages
pages = Any[
    "Home" => "index.md",
    "APIs" => [
        "EmeraldConstants"        => "modules/EmeraldConstants.md",
        "EarthSurface"            => "modules/EarthSurface.md",
        "WaterPhysics"            => "modules/WaterPhysics.md",
        "EmeraldNamespace"        => "modules/EmeraldNamespace.md",
        "LeafOptics"              => "modules/LeafOptics.md",
        "CanopyRadiativeTransfer" => "modules/CanopyRadiativeTransfer.md",
        "Photosynthesis"          => "modules/Photosynthesis.md",
        "SoilHydraulics"          => "modules/SoilHydraulics.md",
        "PlantHydraulics"         => "modules/PlantHydraulics.md",
        "StomataModels"           => "modules/StomataModels.md",
        "SoilPlantAirContinuum"   => "modules/SoilPlantAirContinuum.md"
    ],
    "Tips" => "tips.md",
];

@show pages;


# format the docs
mathengine = MathJax(
    Dict(
        :TeX => Dict(
            :equationNumbers => Dict(:autoNumber => "AMS"),
            :Macros => Dict()
        )
    )
);

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    mathengine = mathengine,
    collapselevel = 1,
    assets = ["assets/favicon.ico"]
);


# build the docs
makedocs(
    sitename = "Land",
    format = format,
    clean = false,
    modules = [Land],
    pages = pages
);


# deploy the docs to Github gh-pages
deploydocs(
    repo = "github.com/CliMA/Land.git",
    target = "build",
    devbranch = "main",
    push_preview = true
);

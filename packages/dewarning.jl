#
# This file is meant to synchronize the Package.jl files to resolve their relative paths.
# To run this file, do
#     julia dewarning.jl
#

pkg_names = [
            "CanopyRadiativeTransfer",
            "EarthSurface",
            "EmeraldConstants",
            "EmeraldNamespace",
            "EmeraldOptics",
            "LeafOptics",
            "Photosynthesis",
            "PlantHydraulics",
            "SoilHydraulics",
            "SoilPlantAirContinuum",
            "StomataModels",
            "WaterPhysics"
];

for pkg_name in pkg_names
    _file_out = "$(@__DIR__)/$(pkg_name).jl/src/$(pkg_name).jl";
    if isfile(_file_out)
        rm(_file_out; force = true);
    end;
end;

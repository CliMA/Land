#
# This file is meant to synchronize the Package.jl files to resolve their relative paths.
# To run this file, do
#     julia remove_warning.jl
#

pkg_names = ["CanopyRadiativeTransfer", "EarthSurface", "EmeraldConstants", "EmeraldNamespace", "LeafOptics", "Photosynthesis", "PlantHydraulics", "SoilHydraulics", "SoilPlantAirContinuum",
             "StomataModels", "WaterPhysics"];

for pkg_name in pkg_names
    _file_in = "$(@__DIR__)/../src/modules/$(pkg_name).jl";
    _file_out = "$(@__DIR__)/$(pkg_name).jl/src/$(pkg_name).jl";

    # synchronize files if both files exist
    @info pkg_name;
    if isfile(_file_in) && isfile(_file_out)
        open(_file_out, "w") do _wfile
            for _line in readlines(_file_in);
                # resolve the relative path
                if occursin("../../packages/$(pkg_name).jl/src/", _line)
                    _line = replace(_line, "../../packages/$(pkg_name).jl/src/" => "");
                end;

                # resolve the interdependencies
                if length(_line) > 8 && _line[1:8] == "using .."
                    _line = "#" * replace(_line, ".." => "");
                end;
                if length(_line) > 9 && _line[1:9] == "import .."
                    _line = "#" * replace(_line, ".." => "");
                end;

                # write line to file
                write(_wfile, _line * "\n");
            end;
        end;
    end;
end;

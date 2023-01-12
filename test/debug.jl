#
#
# This file is meant for debugging and profiling the code. Run the command below to check the memory allocations for the main code
#     julip --track-allocation=user debug.jl
#
#
using Land
using Profile


# create a spac to work on
d_spac = EmeraldNamespace.MonoMLTreeSPAC{Float64}();
d_spac.SOIL.ALBEDO = EmeraldNamespace.BroadbandSoilAlbedo{Float64}();


# clear the memory allocation
Profile.clear_malloc_data();

CanopyRadiativeTransfer.soil_albedo!(d_spac.CANOPY, d_spac.SOIL);


#=
# run the model for 30 times
for _i in 1:30
    SoilPlantAirContinuum.soil_plant_air_continuum!(d_spac, 120.0; p_on = false, t_on = false, θ_on = false);
end;
=#

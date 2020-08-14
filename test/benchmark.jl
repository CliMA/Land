# Benchmark the CanopyRadiation model

using BenchmarkTools, CanopyRadiation, Profile, Revise

FT = Float32;
canopy_rt = create_canopy_rt(FT, nLayer=20);
wl_set    = create_wave_length(FT);
rt_dim    = create_rt_dims(canopy_rt, wl_set);
canRad_rt = create_canopy_rads(FT, rt_dim);
canOpt_rt = create_canopy_opticals(FT, rt_dim);
sunRad_rt = create_incoming_radiation(wl_set);
soil      = create_soil_opticals(wl_set);
angles    = SolarAngles{FT}();
rt_con    = create_rt_container(FT, rt_dim);

arrayOfLeaves = [create_leaf_bios(FT, rt_dim) for i in 1:rt_dim.nLayer];
for i in 1:rt_dim.nLayer
    fluspect!(arrayOfLeaves[i], wl_set);
end

canopy_geometry!(canopy_rt, angles, canOpt_rt, rt_con);
canopy_matrices!(arrayOfLeaves, canOpt_rt);
short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, rt_con);
canopy_fluxes!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, arrayOfLeaves, wl_set, rt_con);
SIF_fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, wl_set, rt_con, rt_dim);

Profile.clear_malloc_data();

#@btime canopy_geometry!($canopy_rt, $angles, $canOpt_rt, $rt_con);
#@btime canopy_matrices!($arrayOfLeaves, $canOpt_rt);
#@btime short_wave!($canopy_rt, $canOpt_rt, $canRad_rt, $sunRad_rt, $soil, $rt_con);
#@btime canopy_fluxes!($canopy_rt, $canOpt_rt, $canRad_rt, $sunRad_rt, $soil, $arrayOfLeaves, $wl_set, $rt_con);
#@btime SIF_fluxes!($arrayOfLeaves, $canOpt_rt, $canRad_rt, $canopy_rt, $soil, $wl_set, $rt_con, $rt_dim);

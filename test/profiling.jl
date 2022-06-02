# Benchmark the CanopyLayers model

using BenchmarkTools, CanopyLayers, Profile

FT      = Float32;
can     = create_canopy_rt(FT, nLayer=20);
wls     = create_wave_length(FT);
rt_dim  = create_rt_dims(can, wls);
can_rad = create_canopy_rads(FT, rt_dim);
can_opt = create_canopy_opticals(FT, rt_dim);
in_rad  = create_incoming_radiation(wls);
soil    = SoilOpticals{FT}(wls);
angles  = SolarAngles{FT}();
rt_con  = create_rt_cache(FT, rt_dim);
leaves  = [create_leaf_bios(FT, rt_dim) for i in 1:rt_dim.nLayer];

canopy_geometry!(can, angles, can_opt, rt_con);
canopy_matrices!(leaves, can_opt);
short_wave!(can, can_opt, can_rad, in_rad, soil, rt_con);
canopy_fluxes!(can, can_opt, can_rad, in_rad, soil, leaves, wls, rt_con);
SIF_fluxes!(leaves, can_opt, can_rad, can, soil, wls, rt_con, rt_dim);

Profile.clear_malloc_data();

#@btime canopy_geometry!($can, $angles, $can_opt, $rt_con);
#@btime canopy_matrices!($leaves, $can_opt);
#@btime short_wave!($can, $can_opt, $can_rad, $in_rad, $soil, $rt_con);
#@btime canopy_fluxes!($can, $can_opt, $can_rad, $in_rad, $soil, $leaves, $wls, $rt_con);
#@btime SIF_fluxes!($leaves, $can_opt, $can_rad, $can, $soil, $wls, $rt_con, $rt_dim);

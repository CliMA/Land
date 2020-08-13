# Benchmark the CanopyRadiation model

using BenchmarkTools, CanopyRadiation, Profile, Revise

FT = Float32;
wl_set    = create_wave_length(FT);
leaf_1    = create_leaf_bios(FT, wl_set.nWL, wl_set.nWLE, wl_set.nWLF);
leaf_2    = create_leaf_bios(FT, wl_set.nWL, wl_set.nWLE, wl_set.nWLF);
canopy_rt = Canopy4RT{FT}(nLayer=20, LAI=FT(3));
canRad_rt = CanopyRads{FT}(nWL=wl_set.nWL, nWLF=wl_set.nWLF, nIncl=length(canopy_rt.litab), nAzi=length(canopy_rt.lazitab), nLayer=canopy_rt.nLayer);
canOpt_rt = create_canopy_opticals(FT, wl_set.nWL, canopy_rt.nLayer, length(canopy_rt.lazitab), length(canopy_rt.litab));
sunRad_rt = create_incoming_radiation(wl_set.sWL);
soil      = create_soil_opticals(wl_set);
angles    = SolarAngles{FT}();
rt_con    = create_rt_container(canopy_rt, canOpt_rt, angles, soil, wl_set);

arrayOfLeaves = [create_leaf_bios(FT, wl_set.nWL, wl_set.nWLE, wl_set.nWLF) for i in 1:canopy_rt.nLayer];
for i in 1:canopy_rt.nLayer
    fluspect!(arrayOfLeaves[i], wl_set);
end

canopy_geometry!(canopy_rt, angles, canOpt_rt, rt_con);
canopy_matrices!(arrayOfLeaves, canOpt_rt);
short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, rt_con);
canopy_fluxes!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, arrayOfLeaves, wl_set, rt_con);
sif_fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, wl_set, rt_con);

Profile.clear_malloc_data();

@btime sif_fluxes!($arrayOfLeaves, $canOpt_rt, $canRad_rt, $canopy_rt, $soil, $wl_set, $rt_con);

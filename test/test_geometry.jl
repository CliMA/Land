# Benchmark the CanopyRadiation model

using BenchmarkTools, CanopyRadiation, Profile, Revise

FT = Float32;
wl_set    = create_wave_length(FT);
leaf_1    = create_leaf_bios(FT, wl_set.nwl, wl_set.nWlE, wl_set.nWlF);
leaf_2    = create_leaf_bios(FT, wl_set.nwl, wl_set.nWlE, wl_set.nWlF);
canopy_rt = Canopy4RT{FT}(nLayer=20, LAI=FT(3));
canRad_rt = CanopyRads{FT}(nWL=wl_set.nwl, nWLf=wl_set.nWlF, nIncl=length(canopy_rt.litab), nAzi=length(canopy_rt.lazitab), nLayer=canopy_rt.nLayer);
canOpt_rt = create_canopy_opticals(FT, wl_set.nwl, canopy_rt.nLayer, length(canopy_rt.lazitab), length(canopy_rt.litab));
sunRad_rt = create_incoming_radiation(wl_set.swl);
soil      = create_soil_opticals(wl_set);
angles    = SolarAngles{FT}();
rt_con    = create_rt_container(canopy_rt, canOpt_rt, angles);

collections = initialize_rt_module(LAI=FT(3));

wl_blue   = FT(450.0);
wl_red    = FT(600.0);
wl_FarRed = FT(740.0);
wl_Red    = FT(685.0);
ind_wle_blue = argmin( abs.(wl_set.wle .- wl_blue  ) );
ind_wle_red  = argmin( abs.(wl_set.wle .- wl_red   ) );
ind_wlf_FR   = argmin( abs.(wl_set.wlf .- wl_FarRed) );
ind_wlf_R    = argmin( abs.(wl_set.wlf .- wl_Red   ) );
ind_red      = argmin( abs.(wl_set.wl  .- wl_Red   ) );
ind_NIR      = argmin( abs.(wl_set.wl  .- 800      ) );

leaf_2.Cab = FT(80.0 );
leaf_2.Cw  = FT(0.012);
fluspect!(leaf_1, wl_set);
fluspect!(leaf_2, wl_set);

arrayOfLeaves = [create_leaf_bios(FT, wl_set.nwl, wl_set.nWlE, wl_set.nWlF) for i in 1:canopy_rt.nLayer];
for i in 1:canopy_rt.nLayer
    fluspect!(arrayOfLeaves[i], wl_set);
end

canopy_geometry!(canopy_rt, angles, canOpt_rt, rt_con);
canopy_matrices!(arrayOfLeaves, canOpt_rt);
short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, rt_con);
canopy_fluxes!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, arrayOfLeaves, wl_set);

Profile.clear_malloc_data();

#@btime short_wave!($canopy_rt, $canOpt_rt, $canRad_rt, $sunRad_rt, $soil, $rt_con);

###############################################################################
#
# Packages notes
# ConstrainedRootSolvers is registered, install it like normal Packages
# CanopyLayers is registered, and will be available 3 days later
# For now install it using
# pkg> add https://github.com/Yujie-W/CanopyLayers.jl
#
###############################################################################
using CSV
using CanopyLayers
using DataFrames
using PyPlot
using Statistics

# initialize the parameters and functions
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
    arrayOfLeaves[i].Cab = 40;
    fluspect!(arrayOfLeaves[i], wl_set);
end


oco3_data = DataFrame!(CSV.File("../data/test_sif.csv"));

# Run code over OCO geometries
SZA = oco3_data[:sza];
RAA = oco3_data[:raa];
VZA = oco3_data[:vza];
LAI = oco3_data[:lai];




reflVIS = [];
reflNIR = [];
SIF_R   = [];
SIF_FR  = [];

for i = 1:length(VZA)
    angles.tts = SZA[i];
    angles.psi = RAA[i];
    angles.tto = VZA[i];

    canopy_rt.LAI  = LAI[i];
    canopy_rt.iLAI = LAI[i] / rt_dim.nLayer;
    canopy_geometry!(canopy_rt, angles, canOpt_rt, rt_con);
    canopy_matrices!(arrayOfLeaves, canOpt_rt);
    short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, rt_con);
    canopy_fluxes!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, arrayOfLeaves, wl_set, rt_con);
    SIF_fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, wl_set, rt_con, rt_dim);

    push!(reflVIS, canRad_rt.alb_obs[28]);
    push!(reflNIR, canRad_rt.alb_obs[52]);
    push!(SIF_R  , (canRad_rt.SIF_obs[19] + canRad_rt.SIF_obs[20]) / 2); # 740 nm, 44 -> 19, not 8
    push!(SIF_FR , canRad_rt.SIF_obs[26]); # 771 nm, 50 -> 26, not 20
end



begin
    figure(1)
    clf();
    subplot(2,2,1);
    plot(reflVIS, oco3_data.rad757, "k+");

    subplot(2,2,2);
    plot(reflNIR, oco3_data.rad771, "k+");

    subplot(2,2,3);
    plot(SIF_R, oco3_data.sif740, "k+");

    subplot(2,2,4);
    plot(SIF_FR, oco3_data.sif771, "k+");
end




@show wl_set.WLF;
@show names(oco3_data);

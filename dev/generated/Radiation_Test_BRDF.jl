# Use Julia Plots package and switch to plotly js option:
using PyPlot

#using Revise
using Parameters

using Land
using Land.CanopyRT

#@unpack FT, leafbio, canopy, angles, canOpt, canRad,sunRad,soil, wl, wle, wlf = CanopyRT

const FT = Float32

wl_set = create_wl_para_set(FT)
leaf = create_leaf_bio(FT, wl_set.nwl, wl_set.nWlE, wl_set.nWlF);
canopy_rt = Canopy4RT{FT, 20, 3.0}()
canRad_rt = CanopyRadiation{FT, wl_set.nwl, wl_set.nWlF, length(canopy_rt.litab), length(canopy_rt.lazitab), canopy_rt.nlayers}()
canOpt_rt = create_canopy_optical(FT, wl_set.nwl, canopy_rt.nlayers, length(canopy_rt.lazitab), length(canopy_rt.litab); using_marray=false)
sunRad_rt = create_incoming_radiation(FT, wl_set.swl);

# show leaf Chlorophyll content:
@show leaf.Cab

# Run Fluspect:
fluspect!(leaf, wl_set);

# Plot Mb matrix
figure()
contourf(wl_set.wle, wl_set.wlf, leaf.Mb)
xlabel("Excitation wavelength (nm)")
ylabel("Emission wavelength (nm)")
title("Fluorescence backward (refl) emission (Mb)")
gcf()

# Plot Mf matrix
figure()
contourf(wl_set.wle, wl_set.wlf, leaf.Mf)
xlabel("Excitation wavelength (nm)")
ylabel("Emission wavelength (nm)")
title("Fluorescence forward (transmission) emission (Mf)")
gcf()

# Define a few wavelengths:
wl_blue = 450.0;
wl_red = 600.0;
wl_FarRed = 740.0;
wl_Red = 685.0;
ind_wle_blue  = argmin(abs.(wl_set.wle .-wl_blue));
ind_wle_red = argmin(abs.(wl_set.wle .-wl_red));
ind_wlf_FR  = argmin(abs.(wl_set.wlf .-wl_FarRed));
ind_wlf_R  = argmin(abs.(wl_set.wlf .-wl_Red));
ind_red = argmin(abs.(wl_set.wl .-wl_Red));
ind_NIR = argmin(abs.(wl_set.wl .-800));

# Plot some cross section in wle and wlf space:
figure()
plot(wl_set.wlf,leaf.Mf[:,ind_wle_blue], color=:black, lw = 2 , label="Forward SIF, excited at $wl_blue nm")
plot(wl_set.wlf,leaf.Mb[:,ind_wle_blue],color=:orange, lw = 2 , label="Backward SIF, excited at $wl_blue nm" )
plot(wl_set.wlf,leaf.Mf[:,ind_wle_red], color=:black,ls="--",lw=2, label="Forward SIF, excited at $wl_red nm" )
plot(wl_set.wlf,leaf.Mb[:,ind_wle_red], color=:orange,ls="--",lw=2, label="Backward SIF, excited at $wl_red nm" )
xlabel("Wavelength (nm)")
ylabel("Fluorescence")
legend()
gcf()

# Plot some cross section in wle and wlf space:
figure()
plot(wl_set.wle,leaf.Mf[ind_wlf_FR,:], color=:black, lw = 2 , label="Forward SIF, emitted at $wl_FarRed nm")
plot(wl_set.wle,leaf.Mb[ind_wlf_FR,:],color=:orange, lw = 2 , label="Backward SIF, emitted at $wl_FarRed nm" )
plot(wl_set.wle,leaf.Mf[ind_wlf_R,:], color=:black,ls="--",lw=2, label="Forward SIF, emitted at $wl_Red nm" )
plot(wl_set.wle,leaf.Mb[ind_wlf_R,:], color=:orange,ls="--",lw=2, label="Backward SIF, emitted at $wl_Red nm" )
xlabel("Absorbed Wavelength (nm)")
ylabel("Fluorescence")
legend()
gcf()

# Let's create a leaf with a different Cab and Cw (water) content (<span style="color:red">Try changing other pigment contents, plot leaf reflectance and transmissions and explain where (spectrally) and why reflectance and transmission changes</span>):
leaf_2 = create_leaf_bio(FT, wl_set.nwl, wl_set.nWlE, wl_set.nWlF);
leaf_2.Cab = 80
leaf_2.Cw = 0.012
# show leaf Chlorophyll content:
# and Run Fluspect:
fluspect!(leaf_2, wl_set);

figure()
plot(wl_set.wl,1 .-leaf.τ_SW, color=:black, lw = 2 , label="Leaf Transmission")
plot(wl_set.wl,leaf.ρ_SW,color=:orange, lw = 2 , label="Leaf Reflectance" )
plot(wl_set.wl,1 .-leaf_2.τ_SW, color=:black, ls="--",lw=2, label="Leaf ##2 Transmission")
plot(wl_set.wl,leaf_2.ρ_SW,color=:orange, ls="--",lw=2, label="Leaf ##2 Reflectance" )
xlabel("Wavelength (nm)")
legend()
gcf()

# This is to be changed later but at the moment, we need to generate an Array of leaves, basically for each layer of the canopy

arrayOfLeaves = [create_leaf_bio(FT, wl_set.nwl, wl_set.nWlE, wl_set.nWlF) for i in 1:canopy_rt.nlayers]
for i in 1:canopy_rt.nlayers
    fluspect!(arrayOfLeaves[i],  wl_set)
end

# Set Soil albedo to 0.2
soil = SoilOpti{FT}(wl_set.wl, FT(0.2)*ones(FT, length(wl_set.wl)), FT[0.1], FT(290.0))
angles = SolarAngles{FT}()

# Compute Canopyoptical properties dependend on sun-sensor and leaf angle distributions:
compute_canopy_geometry!(canopy_rt, angles, canOpt_rt)
# Compute RT matrices with leaf reflectance and transmissions folded in:
compute_canopy_matrices!(arrayOfLeaves, canOpt_rt);
# Perform SW radiation transfer:
simulate_short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil);
# Compute outgoing SIF flux (using constant fluorescence efficiency at the chloroplast level)
derive_canopy_fluxes!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, arrayOfLeaves, wl_set);

SIF_FR = Float32[]
SIF_R = Float32[]
reflVIS = Float32[]
reflNIR = Float32[]

# Just running the code over all geometries:

# Set sun SZA to 30 degrees
angles.tts=30
# Set 0 azimuth (principal plane)
angles.psi=0
# LAI of 3:
canopy_rt.LAI = 3
# Define VZA
VZA=collect(-89.5:0.5:89.5)

for VZA_ in VZA
    angles.tto=VZA_
    compute_canopy_geometry!(canopy_rt, angles, canOpt_rt)
    compute_canopy_matrices!(arrayOfLeaves, canOpt_rt);
    simulate_short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil);
    computeSIF_Fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, wl_set);
    # Handpicked indices in
    push!(reflVIS, canRad_rt.alb_obs[ind_red])
    push!(reflNIR, canRad_rt.alb_obs[ind_NIR])
    push!(SIF_R , canRad_rt.SIF_obs[ind_wlf_R])
    push!(SIF_FR, canRad_rt.SIF_obs[ind_wlf_FR ])
end

# Plots Visible
figure()
plot(VZA, reflVIS, label="Red Reflectance", lw=2)
plot(VZA, SIF_R/30, label="Red SIF (/30)", lw=2)
xlabel("Viewing Zenith Angle")
legend()
gcf()

# Plot Near Infrared:
figure()
plot(VZA, reflNIR, label="NIR Reflectance", lw=2)
plot(VZA, SIF_FR/6, label="Far Red SIF (/6)", lw=2)
xlabel("Viewing Zenith Angle")
gcf()

reflVIS = Float32[]
reflNIR = Float32[]
SIF_FR = Float32[]
SIF_R  = Float32[]
angles.tts=48
angles.psi=0
canopy_rt.LAI=3.22
for psi=0:360
    angles.psi=psi
    for VZA=0:1:85
        angles.tto=VZA

        compute_canopy_geometry!(canopy_rt, angles, canOpt_rt)
        compute_canopy_matrices!(arrayOfLeaves, canOpt_rt);
        simulate_short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil);
        computeSIF_Fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, wl_set);
        push!(reflVIS, canRad_rt.alb_obs[28])
        push!(reflNIR, canRad_rt.alb_obs[52])
        push!(SIF_R , canRad_rt.SIF_obs[8])
        push!(SIF_FR, canRad_rt.SIF_obs[20])
    end
end

A = reshape(reflNIR, ( 86,361));
B = reshape(reflVIS, ( 86,361));
SIFFER = reshape(SIF_R, ( 86,361));
SIFFER_FR = reshape(SIF_FR, ( 86,361));

##heatmap(A, cmap=)
figure()
subplot(1,1,1, polar=true)
grid(false)
hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  A,  cmap=:viridis)
title("NIR reflectance BRDF")
colorbar()
gcf()

##heatmap(A, cmap=)
figure()
subplot(1,1,1, polar=true)
grid(false)
hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  B,  cmap=:viridis)
title("VIS reflectance BRDF")
colorbar()
gcf()

figure()
subplot(1,1,1, polar=true)
grid(false)
hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  SIFFER, cmap=:viridis)
title("Red SIF emission BRDF")
colorbar()
gcf()

figure()
subplot(1,1,1, polar=true)
grid(false)
hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  SIFFER_FR, cmap=:viridis)
title("Far Red SIF emission BRDF")
colorbar()
gcf()

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl


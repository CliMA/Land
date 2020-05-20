# Use Julia Plots package and switch to plotly js option:
using Plots
#gr()
#using Conda
#Conda.add("matplotlib")

pyplot()
#plotlyjs()

#using Revise
using Parameters

using Land
using Land.CanopyRT

#@unpack FT, leafbio, canopy, angles, canOpt, canRad,sunRad,soil, wl, wle, wlf = CanopyRT

leaf = leafbio{FT, length(wl), length(wle), length(wlf),length(wle)*length(wlf)}();
# show leaf Chlorophyll content:
@show leaf.Cab

# Run Fluspect:
CanopyRT.fluspect!(leaf, CanopyRT.optis);

# Plot Mb matrix
contourf(wle, wlf, leaf.Mb)
xlabel!("Excitation wavelength (nm)")
ylabel!("Emission wavelength (nm)")
title!("Fluorescence backward (refl) emission (Mb)")

# Plot Mf matrix
contourf(wle, wlf, leaf.Mf)
xlabel!("Excitation wavelength (nm)")
ylabel!("Emission wavelength (nm)")
title!("Fluorescence forward (transmission) emission (Mf)")

# Define a few wavelengths:
wl_blue = 450.0;
wl_red = 600.0;
wl_FarRed = 740.0;
wl_Red = 685.0;
ind_wle_blue  = argmin(abs.(wle .-wl_blue));
ind_wle_red = argmin(abs.(wle .-wl_red));
ind_wlf_FR  = argmin(abs.(wlf .-wl_FarRed));
ind_wlf_R  = argmin(abs.(wlf .-wl_Red));
ind_red = argmin(abs.(wl .-wl_Red));
ind_NIR = argmin(abs.(wl .-800));

# Plot some cross section in wle and wlf space:
plot(wlf,leaf.Mf[:,ind_wle_blue], color=[:black], lw = 2 , label="Forward SIF, excited at $wl_blue nm")
plot!(wlf,leaf.Mb[:,ind_wle_blue],color=[:orange], lw = 2 , label="Backward SIF, excited at $wl_blue nm" )
plot!(wlf,leaf.Mf[:,ind_wle_red], color=[:black],line=(:dash,2), label="Forward SIF, excited at $wl_red nm" )
plot!(wlf,leaf.Mb[:,ind_wle_red], color=[:orange],line=(:dash,2), label="Backward SIF, excited at $wl_red nm" )
xlabel!("Wavelength (nm)")
ylabel!("Fluorescence")

# Plot some cross section in wle and wlf space:
plot(wle,leaf.Mf[ind_wlf_FR,:], color=[:black], lw = 2 , label="Forward SIF, emitted at $wl_FarRed nm")
plot!(wle,leaf.Mb[ind_wlf_FR,:],color=[:orange], lw = 2 , label="Backward SIF, emitted at $wl_FarRed nm" )
plot!(wle,leaf.Mf[ind_wlf_R,:], color=[:black],line=(:dash,2), label="Forward SIF, emitted at $wl_Red nm" )
plot!(wle,leaf.Mb[ind_wlf_R,:], color=[:orange],line=(:dash,2), label="Backward SIF, emitted at $wl_Red nm" )
xlabel!("Absorbed Wavelength (nm)")
ylabel!("Fluorescence")

# Let's create a leaf with a different Cab and Cw (water) content (<span style="color:red">Try changing other pigment contents, plot leaf reflectance and transmissions and explain where (spectrally) and why reflectance and transmission changes</span>):
leaf_2 = leafbio{FT, length(wl), length(wle), length(wlf),length(wle)*length(wlf)}();
leaf_2.Cab = 80
leaf_2.Cw = 0.012
# show leaf Chlorophyll content:
# and Run Fluspect:
CanopyRT.fluspect!(leaf_2, CanopyRT.optis);

plot(wl,1 .-leaf.τ_SW, color=[:black], lw = 2 , label="Leaf Transmission")
plot!(wl,leaf.ρ_SW,color=[:orange], lw = 2 , label="Leaf Reflectance" )
plot!(wl,1 .-leaf_2.τ_SW, color=[:black], line=(:dash,2), label="Leaf ##2 Transmission")
plot!(wl,leaf_2.ρ_SW,color=[:orange], line=(:dash,2), label="Leaf ##2 Reflectance" )
xlabel!("Wavelength (nm)")

# This is to be changed later but at the moment, we need to generate an Array of leaves, basically for each layer of the canopy
arrayOfLeaves = Array{leafbio{FT,length(wl), length(wle), length(wlf),length(wle)*length(wlf)}, 1}(undef, CanopyRT.canopy.nlayers)
for i = 1:CanopyRT.canopy.nlayers
    ##@show i
    arrayOfLeaves[i] = leafbio{FT, length(wl), length(wle), length(wlf),length(wle)*length(wlf)}()
    CanopyRT.fluspect!(arrayOfLeaves[i], CanopyRT.optis)
end

# Set Soil albedo to 0.2
CanopyRT.soil.albedo_SW[:] .=0.2;
# Compute Canopyoptical properties dependend on sun-sensor and leaf angle distributions:
CanopyRT.computeCanopyGeomProps!(canopy, angles,canOpt);
# Compute RT matrices with leaf reflectance and transmissions folded in:
CanopyRT.computeCanopyMatrices!(arrayOfLeaves,canOpt);
# Perform SW radiation transfer:
CanopyRT.RTM_SW!(canopy, canOpt, canRad,sunRad, CanopyRT.soil);
# Compute outgoing SIF flux (using constant fluorescence efficiency at the chloroplast level)
CanopyRT.computeSIF_Fluxes!(arrayOfLeaves, canOpt, canRad, canopy, CanopyRT.soil);

SIF_FR = Float32[]
SIF_R = Float32[]
reflVIS = Float32[]
reflNIR = Float32[]

# Just running the code over all geometries:

# Set sun SZA to 30 degrees
CanopyRT.angles.tts=30
# Set 0 azimuth (principal plane)
CanopyRT.angles.psi=0
# LAI of 3:
CanopyRT.canopy.LAI = 3
# Define VZA
VZA=collect(-89.5:0.5:89.5)

for VZA_ in VZA
    CanopyRT.angles.tto=VZA_
    CanopyRT.computeCanopyGeomProps!(canopy, angles,canOpt);
    CanopyRT.computeCanopyMatrices!(arrayOfLeaves,canOpt);
    CanopyRT.RTM_SW!(canopy, canOpt, canRad,sunRad, CanopyRT.soil);
    CanopyRT.computeSIF_Fluxes!(arrayOfLeaves, canOpt, canRad, canopy, CanopyRT.soil);
    # Handpicked indices in
    push!(reflVIS, canRad.alb_obs[ind_red])
    push!(reflNIR, canRad.alb_obs[ind_NIR])
    push!(SIF_R , canRad.SIF_obs[ind_wlf_R])
    push!(SIF_FR, canRad.SIF_obs[ind_wlf_FR ])
end

# Plots Visible
plot(VZA, reflVIS, label="Red Reflectance", lw=2)
plot!(VZA, SIF_R/30, label="Red SIF (/30)", lw=2)
xlabel!("Viewing Zenith Angle")

# Plot Near Infrared:
plot(VZA, reflNIR, label="NIR Reflectance", lw=2)
plot!(VZA, SIF_FR/6, label="Far Red SIF (/6)", lw=2)
xlabel!("Viewing Zenith Angle")

reflVIS = Float32[]
reflNIR = Float32[]
SIF_FR = Float32[]
SIF_R  = Float32[]
CanopyRT.angles.tts=48
CanopyRT.angles.psi=0
CanopyRT.canopy.LAI=3.22
for psi=0:360
    CanopyRT.angles.psi=psi
    for VZA=0:1:85
        CanopyRT.angles.tto=VZA

        CanopyRT.computeCanopyGeomProps!(canopy, angles,canOpt);
        CanopyRT.computeCanopyMatrices!(arrayOfLeaves,canOpt);
        CanopyRT.RTM_SW!(canopy, canOpt, canRad,sunRad, CanopyRT.soil);
        CanopyRT.computeSIF_Fluxes!(arrayOfLeaves, canOpt, canRad, canopy, CanopyRT.soil);
        push!(reflVIS, canRad.alb_obs[28])
        push!(reflNIR, canRad.alb_obs[52])
        push!(SIF_R , canRad.SIF_obs[8])
        push!(SIF_FR, canRad.SIF_obs[20])
    end
end

A = reshape(reflNIR, ( 86,361));
B = reshape(reflVIS, ( 86,361));
SIFFER = reshape(SIF_R, ( 86,361));
SIFFER_FR = reshape(SIF_FR, ( 86,361));

##heatmap(A, cmap=)
hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  A,  proj=:polar, color=:viridis, alpha=0.5)
title!("NIR reflectance BRDF")

##heatmap(A, cmap=)
hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  B,  proj=:polar, color=:viridis)
title!("VIS reflectance BRDF")

hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  SIFFER, proj=:polar, color=:viridis)
title!("Red SIF emission BRDF")

hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  SIFFER_FR, proj=:polar, color=:viridis)
title!("Far Red SIF emission BRDF")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl


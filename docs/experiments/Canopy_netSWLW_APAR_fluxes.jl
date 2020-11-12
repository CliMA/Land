
# Add PATH (adapt!)
push!(LOAD_PATH, "/Volumes/cfranken/code/gitHub/LSM-SPAM/src/Leaf/");
push!(LOAD_PATH, "/Volumes/cfranken/code/gitHub/LSM-SPAM/src/Utils/");
#push!(LOAD_PATH, "/Volumes/cfranken/code/gitHub/LSM-SPAM/src/tools/");

# Use Plots:
using Plots
using BenchmarkTools
using Statistics
using Parameters

using Revise
using CanopyRTMod
using LeafPhotosynthesis

@unpack wl,wle,wlf, soil = CanopyRTMod;
# Create an array of standard leaves (needs to be in Module later on:
arrayOfLeaves = Array{leafbio{FT,length(wl), length(wle), length(wlf),length(wle)*length(wlf)}, 1}(undef, CanopyRTMod.canopy.nlayers)
for i = 1:CanopyRTMod.canopy.nlayers
    #@show i
    arrayOfLeaves[i] = leafbio{FT, length(wl), length(wle), length(wlf),length(wle)*length(wlf)}()
    CanopyRTMod.fluspect!(arrayOfLeaves[i], CanopyRTMod.optis)
end


# 4 Different steps to compute Short-Wave RT
@time CanopyRTMod.computeCanopyGeomProps!(canopy, angles,canOpt)
@time CanopyRTMod.computeCanopyMatrices!(arrayOfLeaves,canOpt);
@time CanopyRTMod.RTM_SW!(canopy, canOpt, canRad,sunRad, CanopyRTMod.soil);
@time CanopyRTMod.deriveCanopyFluxes!(canopy, canOpt, canRad,sunRad, CanopyRTMod.soil, arrayOfLeaves);
# Compute Long Wave (Last term is LW incoming in W m^-2)
@time CanopyRTMod.computeThermalFluxes!(arrayOfLeaves, canOpt, canRad, canopy, soil, [Float32(400.0)]);

#@show arrayOfLeaves[1].kChlrel
# Layer Temperatures are here:
@show canRad.T_sun;
@show canRad.T_shade;

# Net Energy fluxes
@show canRad.intNetLW_shade;
@show canRad.intNetLW_sunlit;
@show canRad.intNetSW_shade;
@show canRad.intNetSW_sunlit;
@show canRad.RnSoilLW;
@show canRad.RnSoil;
#@show 1e6 * canRad.absPAR_sunCab;
@show 1e6 * canRad.absPAR_shadeCab;

plot(wl, canRad.netSW_sunlit)

iLAI = canopy.LAI/canopy.nlayers
plot(wl,1 ./canOpt.Ps[10] / iLAI * canRad.netSW_sunlit[:,10])
plot!(wl, sunRad.E_direct)

using Leaf
l = leaf_params{Float32}();
l2 = leaf_params{Float32}();
# Create a Flux structure
f = LeafPhotosynthesis.fluxes{Float32}();
l.vcmax25=120
l.jmax25=l.vcmax25*1.8

# initialize some reasonable values
f.je = 100;
f.gbc = 100;
f.gbv = 100;
f.ceair=1500;
f.eair = 1500;
f.APAR = 1.0e6 * canRad.absPAR_shadeCab[1];
1e6*mean(canRad.absPAR_sunCab[:,:,1])

A = similar(canRad.absPAR_sunCab)
Ashade = similar(canRad.absPAR_shadeCab)
I = CartesianIndices(A)
IShade = CartesianIndices(Ashade)
for i in I
    LeafPhotosynthesis.LeafPhotosynthesis(f,l,Float32(298.0), Float32(1.0e6) * canRad.absPAR_sunCab[i]);
    A[i]=f.an;
end
for i in IShade
    LeafPhotosynthesis.LeafPhotosynthesis(f,l,Float32(298.0), Float32(1.0e6) * canRad.absPAR_shadeCab[i]);
    Ashade[i]=f.an;
end

LeafPhotosynthesis.LeafPhotosynthesis(f,l,Float32(298.0), Float32(1.0e6) * mean(canRad.absPAR_sunCab[:,:,1]));

using Plots

plot(Float32(1.0e6) * canRad.absPAR_sunCab[:,:,1][:], A[:,:,1][:], seriestype = :scatter,legend=:bottomright, label="Individual sunlit leaves in one layer  (A(APAR))")
plot!([Float32(1.0e6) * mean(canRad.absPAR_sunCab[:,:,1])], [f.an], seriestype = :scatter,  label="A(average(APAR))")
plot!([Float32(1.0e6) * mean(canRad.absPAR_sunCab[:,:,1])], [mean(A)], seriestype = :scatter, label="average(A((APAR)))")
plot!(Float32(1.0e6) * canRad.absPAR_shadeCab[:], Ashade, seriestype = :scatter,legend=:bottomright, label="Individual shaded layers")



f



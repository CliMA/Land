
# Add PATH (adapt!)
push!(LOAD_PATH, "/Volumes/cfranken/code/gitHub/LSM-SPAM/src/Leaf/");
push!(LOAD_PATH, "/Volumes/cfranken/code/gitHub/LSM-SPAM/src/Utils/");
#push!(LOAD_PATH, "/Volumes/cfranken/code/gitHub/LSM-SPAM/src/tools/");

# Use Plots:
using Plots
#using PhotoStructs
pyplot()

using CanopyRTMod
using BenchmarkTools
using Statistics

wl = CanopyRTMod.wl;
wle = CanopyRTMod.wle;
wlf = CanopyRTMod.wlf;

arrayOfLeaves = Array{leafbio{FT,length(wl), length(wle), length(wlf),length(wle)*length(wlf)}, 1}(undef, CanopyRTMod.canopy.nlayers)
for i = 1:CanopyRTMod.canopy.nlayers
    #@show i
    arrayOfLeaves[i] = leafbio{FT, length(wl), length(wle), length(wlf),length(wle)*length(wlf)}()
    CanopyRTMod.fluspect!(arrayOfLeaves[i], CanopyRTMod.optis)
end
# This actually almost takes 1ms!
#CanopyRTMod.fluspect!(leaf, CanopyRTMod.optis)


CanopyRTMod.computeCanopyGeomProps!(canopy, angles,canOpt)

CanopyRTMod.computeCanopyMatrices!(arrayOfLeaves,canOpt);

@btime CanopyRTMod.RTM_SW!(canopy, canOpt, canRad,sunRad, CanopyRTMod.soil);

@time CanopyRTMod.computeSIF_Fluxes!(arrayOfLeaves, canOpt, canRad, canopy, CanopyRTMod.soil);
#piLo1, piLo2, piLo3, piLo4, F⁻,F⁺,S⁻,S⁺,piLs, piLd =
#println(size(piLs))
#plot(wlf, piLs)
plot(wlf, canRad.SIF_hemi/pi, label="Hemispheric Fluorescence (/pi)")
plot!(wlf, canRad.SIF_obs_sunlit, label="Sunlit direct")
plot!(wlf, canRad.SIF_obs_shaded, label="Shaded direct")
plot!(wlf, canRad.SIF_obs_scattered, label="Scattered")
plot!(wlf, canRad.SIF_obs_soil, label="Soil")

plot(wlf, canRad.SIF_hemi, label="Outgoing hemispheric")
plot!(wlf, canRad.SIF_sum, label="Sum of all layer SIF sources")


plot(wlf, canRad.SIF_hemi./canRad.SIF_sum, label="Outgoing hemispheric/sources")
plot!(wlf, pi*canRad.SIF_obs./canRad.SIF_sum, label="π*SIF_obs/sources")
#plot!(wlf, , label="Sum of all layer SIF sources")

plot(wlf, arrayOfLeaves[1].Mb[:,1])

iLAI = canopy.LAI/canopy.nlayers
ϵ = zeros(20).+0.98
fSun = (canOpt.Ps[1:20]+canOpt.Ps[2:21])/2
wlii = [10]
S⁺,S⁻ = CanopyRTMod.computeThermalFluxes(canRad.T_shade, canRad.T_sun, ϵ, iLAI, canopy.lidf, fSun, wlii)
#plot(S⁺)

@btime CanopyRTMod.deriveCanopyFluxes!(canopy, canOpt, canRad,sunRad, CanopyRTMod.soil, arrayOfLeaves);

plot(wl,canRad.alb_direct, label="Direct hemispheric albedo")
plot!(wl,canRad.alb_diffuse, label="Diffuse hemispheric albedo")
plot!(wl,canRad.alb_obs, label="observed directional (nadir) albedo")


sumLAI = [0:iLAI:canopy.LAI;]
plot(canOpt.Pso,-sumLAI, label="Pso")
plot!(canOpt.Ps,-sumLAI, label="Ps")
plot!(canOpt.Po,-sumLAI, label="Po")
ylabel!("-Cumulative LAI")
xlabel!("Probability ")
title!("solar Ps, outgoing Po, in-out Pso likelihoods")

solar_in = sunRad.E_diffuse+sunRad.E_direct;
soil_absorbed = canRad.E_down[:,end]+canOpt.Es_[:,end]-canRad.E_up[:,end]
canopy_absorbed_diff = sum(canRad.netSW_shade, dims=2)[:,1]
canopy_absorbed_dir = sum(canRad.netSW_sunlit, dims=2)[:,1]
plot(wl,solar_in-CanopyRTMod.canRad.E_up[:,1].+10, label="Net TOA incoming (+10)")
plot!(wl,soil_absorbed, label="Absorbed  by soil")
plot!(wl,soil_absorbed+canopy_absorbed_diff, label="Absorbed  by soil + canopy (diffuse)")
plot!(wl,soil_absorbed+canopy_absorbed_diff+canopy_absorbed_dir, label="Absorbed  by soil + canopy (total)")
xlabel!("Wavelength (nm)")
ylabel!("Radiance (mW m-2 μm-1)")

sum(soil_absorbed+canopy_absorbed_diff+canopy_absorbed_dir-solar_in+CanopyRTMod.canRad.E_up[:,1], dims=1)

plot(wl, arrayOfLeaves[10].ρ_SW, label="Reflectance")
plot!(wl, 1 .-arrayOfLeaves[10].τ_SW, label="1-Transmission")

println("Net Soil radiation (direct): ", canRad.RnSoil_direct, "W/m2")
println("Net Soil radiation (diffuse): ", canRad.RnSoil_diffuse, "W/m2")

println("Incoming direct PAR: ", 1e6*canRad.incomingPAR_direct, "moles m^-2 s^-1")
println("Incoming diffuse PAR: ", 1e6*canRad.incomingPAR_diffuse, "moles m^-2 s^-1")

nl = CanopyRTMod.canopy.nlayers
#normi = 1/mean(canOpt.fs'*canopy.lidf)
plot(sumLAI[1:nl], 1e6*canRad.absPAR_shadeCab, label="Diffuse absorbed PAR")
plot!(sumLAI[1:nl], 1e6*canRad.absPAR_sunCab[5,4,:].*canOpt.Ps[1:nl], label="Direct absorbed PAR*Ps")
plot!(sumLAI[1:nl], 1e6*canRad.absPAR_sunCab[5,4,:], label="Direct absorbed PAR")

mean(canOpt.Ps)

contourf(1e6*canRad.absPAR_sunCab[:,:,1])
#plot!(1e6*canRad.absPAR_sun[:,:,1][:])

# plot(wl[CanopyRTMod.iPAR], arrayOfLeaves[1].kChlrel[CanopyRTMod.iPAR], label="Relative absorbed light by Cab+Car")
plot(wl[CanopyRTMod.iPAR], arrayOfLeaves[1].kChlrel_old[CanopyRTMod.iPAR], label="Relative absorbed light by Cab")


nwl2,nl2 = size(canOpt.R_dd)

nwl2



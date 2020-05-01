
# Add PATH (adapt!)
push!(LOAD_PATH, "/Volumes/cfranken/code/gitHub/Land/src/");
#push!(LOAD_PATH, "/Volumes/cfranken/code/gitHub/Land/src/Utils/");
#push!(LOAD_PATH, "/Volumes/cfranken/code/gitHub/LSM-SPAM/src/tools/");

# Use Plots:
using Plots
#using PhotoStructs
pyplot()

using Revise
using Land
using BenchmarkTools
using Statistics
using Parameters

@unpack CanopyRTMod = Land
@unpack FT,leafbio, canopy, angles, canOpt, canRad,sunRad,soil = Land.CanopyRTMod

wl = CanopyRTMod.wl;
wle = CanopyRTMod.wle;
wlf = CanopyRTMod.wlf;


#@time CanopyRTMod.optis.nr[3:20] .= 2

leaf = leafbio{FT, length(wl), length(wle), length(wlf),length(wle)*length(wlf)}();


using Parameters
@unpack N, Cab, Car, Ant, Cs, Cw, Cm, ρ_SW, τ_SW = CanopyRTMod.leaf;

@btime CanopyRTMod.fluspect!(CanopyRTMod.leaf, CanopyRTMod.optis);

contourf(wle, wlf, CanopyRTMod.leaf.Mb)

contourf(wle, wlf, CanopyRTMod.leaf.Mf)

arrayOfLeaves = Array{leafbio{FT,length(wl), length(wle), length(wlf),length(wle)*length(wlf)}, 1}(undef, CanopyRTMod.canopy.nlayers)
for i = 1:CanopyRTMod.canopy.nlayers
    #@show i
    arrayOfLeaves[i] = leafbio{FT, length(wl), length(wle), length(wlf),length(wle)*length(wlf)}()
    CanopyRTMod.fluspect!(arrayOfLeaves[i], CanopyRTMod.optis)
end

CanopyRTMod.angles.tto=-88
CanopyRTMod.angles.psi=50
@time CanopyRTMod.computeCanopyGeomProps!(canopy, angles,canOpt);
@time CanopyRTMod.computeCanopyMatrices!(arrayOfLeaves,canOpt);
@time CanopyRTMod.RTM_SW!(canopy, canOpt, canRad,sunRad, CanopyRTMod.soil);
CanOpt2 = deepcopy(canOpt)
@time CanopyRTMod.computeSIF_Fluxes!(arrayOfLeaves, canOpt, canRad, canopy, CanopyRTMod.soil);

# Just a test
@time CanopyRTMod.compCanopyOptsExact!(arrayOfLeaves,canOpt, canopy.lidf);

plot(CanOpt2.w[:,10])
plot!(CanOpt2.vb[:,1])
plot!(canOpt.w[:,1]/36)
plot!(canOpt.vb[:,1]/36)

CanopyRTMod.soil.albedo_SW[:] .=0.2;

reflVIS = []
reflNIR = []
SIF_FR = []
SIF_FR1 = []
SIF_FR2 = []
SIF_FR3 = []
SIF_FR4 = []
SIF_R1 = []
SIF_R2 = []
SIF_R3 = []
SIF_R4 = []
SIF_R  = []
Pso = []
Po = []
Ps = []
ko = []
CanopyRTMod.angles.tts=30
CanopyRTMod.angles.psi=0
CanopyRTMod.canopy.LAI = 3
for VZA=-89.5:0.5:89.5
    CanopyRTMod.angles.tto=VZA
    CanopyRTMod.computeCanopyGeomProps!(canopy, angles,canOpt);
    CanopyRTMod.computeCanopyMatrices!(arrayOfLeaves,canOpt);
    CanopyRTMod.RTM_SW!(canopy, canOpt, canRad,sunRad, CanopyRTMod.soil);
    CanopyRTMod.computeSIF_Fluxes!(arrayOfLeaves, canOpt, canRad, canopy, CanopyRTMod.soil);
    push!(reflVIS, canRad.alb_obs[28])
    push!(reflNIR, canRad.alb_obs[52])
    push!(SIF_R , canRad.SIF_obs[8])
    push!(SIF_FR, canRad.SIF_obs[20])
    push!(SIF_FR1, canRad.SIF_obs_sunlit[20])
    push!(SIF_FR2, canRad.SIF_obs_scattered[20])
    push!(SIF_FR3, canRad.SIF_obs_soil[20])
    push!(SIF_FR4, canRad.SIF_obs_shaded[20])
    push!(SIF_R1, canRad.SIF_obs_sunlit[8])
    push!(SIF_R2, canRad.SIF_obs_scattered[8])
    push!(SIF_R3, canRad.SIF_obs_soil[8])
    push!(SIF_R4, canRad.SIF_obs_shaded[8])
    push!(Pso, canOpt.Pso[1])
    push!(Ps, canOpt.Ps[1])
    push!(Po, canOpt.Po[1])
    push!(ko, canOpt.ko[1])
end

CanopyRTMod.angles.tto=89
CanopyRTMod.computeCanopyGeomProps!(canopy, angles,canOpt);
CanopyRTMod.computeCanopyMatrices!(arrayOfLeaves,canOpt);
CanopyRTMod.RTM_SW!(canopy, canOpt, canRad,sunRad, CanopyRTMod.soil);
CanopyRTMod.computeSIF_Fluxes!(arrayOfLeaves, canOpt, canRad, canopy, CanopyRTMod.soil);
plot(canRad.SIF_hemi)
plot!(canRad.SIF_sum)

VZA=-89.5:0.5:89.5
plot(VZA, ko, label="ko")

@show wl[28]
@show wl[52]
CanopyRTMod.angles.tts = 48

plot(VZA, Pso, label="Pso")
plot!(VZA, Po, label="Po")
plot!(VZA, Ps, label="Ps")
#@show wlf


plot(VZA, reflVIS, label="Red Reflectance")
plot!(VZA, SIF_FR/100, label="Red SIF")

plot(VZA, reflNIR, label="NIR Reflectance")
plot!(VZA, SIF_FR, label="Far Red SIF")

plot(VZA, SIF_FR1, label="FR SIF Sunlit")
plot!(VZA, SIF_FR2, label="FR SIF scattered")
plot!(VZA, SIF_FR3, label="FR SIF soil")
plot!(VZA, SIF_FR4, label="FR SIF shaded")
plot!(VZA, SIF_FR, label="FR SIF total")
ylims!(0,5)

plot(VZA, SIF_R1, label="R SIF Sunlit")
plot!(VZA, SIF_R2, label="R SIF scattered")
plot!(VZA, SIF_R3, label="R SIF soil")
plot!(VZA, SIF_R4, label="R SIF shaded")
plot!(VZA, SIF_R, label="R SIF ")
ylims!(0,2)

plot(VZA, SIF_FR./reflNIR./mean(SIF_FR./reflNIR), label="FR/refl")
plot!(VZA, SIF_R./reflVIS./mean(SIF_R./reflVIS), label="R/refl")
#plot!(VZA, , label="Far Red SIF")

CanopyRTMod.soil.albedo_SW[:] .=0.2;
reflRed_SZA = []
reflNIR_SZA = []
CanopyRTMod.canopy.Ω = 1.0
CanopyRTMod.angles.tto=0.2
CanopyRTMod.canopy.LAI=2.
for SZA=0.0:1:75
    CanopyRTMod.angles.tts=SZA
    CanopyRTMod.computeCanopyGeomProps!(canopy, angles,canOpt);
    CanopyRTMod.computeCanopyMatrices!(arrayOfLeaves,canOpt);
    CanopyRTMod.RTM_SW!(canopy, canOpt, canRad,sunRad, CanopyRTMod.soil);
    push!(reflRed_SZA, canRad.alb_direct[28])
    push!(reflNIR_SZA, canRad.alb_direct[52])
end

SZA=0:1:75
plot(SZA, reflNIR_SZA)

SZA=0:1:75
plot(SZA, reflRed_SZA)

# Test plots from Christiaan's papers

reflVIS = Float32[]
reflNIR = Float32[]
SIF_FR = Float32[]
SIF_R  = Float32[]
CanopyRTMod.angles.tts=48
CanopyRTMod.angles.psi=0
CanopyRTMod.canopy.LAI=3.22
for psi=0:360
    CanopyRTMod.angles.psi=psi
    for VZA=0:1:85
        CanopyRTMod.angles.tto=VZA

        CanopyRTMod.computeCanopyGeomProps!(canopy, angles,canOpt);
        CanopyRTMod.computeCanopyMatrices!(arrayOfLeaves,canOpt);
        CanopyRTMod.RTM_SW!(canopy, canOpt, canRad,sunRad, CanopyRTMod.soil);
        CanopyRTMod.computeSIF_Fluxes!(arrayOfLeaves, canOpt, canRad, canopy, CanopyRTMod.soil);
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

using Plots; pyplot()

#heatmap(A, cmap=)
hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  A,  proj=:polar, color=:viridis, alpha=0.5)

#heatmap(A, cmap=)
hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  B,  proj=:polar, color=:viridis)

hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  SIFFER, proj=:polar, color=:viridis)

hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  SIFFER_FR, proj=:polar, color=:viridis)

plot(CanopyRTMod.litab, cumsum(canopy.lidf))

litab_bnd  = FT[[0.,10.,20.,30.,40.,50.,60.,70.,80.,82.,84.,86.,88.] [10.,20.,30.,40.,50.,60.,70.,80.,82.,84.,86.,88., 90.] ];

cumsum(canopy.lidf)



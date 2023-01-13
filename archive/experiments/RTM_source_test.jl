#=
# Add PATH (adapt!)
push!(LOAD_PATH, "/Volumes/cfranken/code/gitHub/LSM-SPAM/src/Leaf/");
push!(LOAD_PATH, "/Volumes/cfranken/code/gitHub/LSM-SPAM/src/Utils/");
#push!(LOAD_PATH, "/Volumes/cfranken/code/gitHub/LSM-SPAM/src/tools/");

# Use Plots:
using Plots
#using PhotoStructs
#plotly()

using CanopyRTMod
using BenchmarkTools
using Statistics
using Parameters
using Revise

@unpack wl,wle,wlf, soil = CanopyRTMod;

leaf = leafbio{FT, length(wl), length(wle), length(wlf)}();
soilAlbedo = 0.07;
incomingLW = 400.0 #W/m2
arrayOfLeaves = Array{leafbio{FT}, 1}(undef, 1)
arrayOfLeaves[1]=leaf

leaf.fqe=0.01
CanopyRTMod.fluspect!(leaf, CanopyRTMod.optis)
CanopyRTMod.computeCanopyGeomProps!(canopy, angles,canOpt)
CanopyRTMod.computeCanopyMatrices!(arrayOfLeaves,canOpt);

using PhysCon
consts = PhysCon.phys{Float32}()

# Hack Thermal optical properties per layer here:
nl = 20
nwl = 1
iLAI    = 0.2;#canopy.LAI/nl;
sigf = canOpt.ddf*leaf.ρ_LW + canOpt.ddb*leaf.τ_LW
sigb = canOpt.ddb*leaf.ρ_LW + canOpt.ddf*leaf.τ_LW
τ_dd = (1 - (1-sigf)*iLAI)*ones(nwl,nl)
ρ_dd = (sigb*iLAI)*ones(nwl,nl)

1.0-τ_dd[1,1]-ρ_dd[1,1]

S⁺ = zeros(nwl,nl).+152*iLAI*0.98
S⁻ = zeros(nwl,nl).+152*iLAI*0.98
S⁺[1,9:11].=200*iLAI
S⁻[1,9:11].=200*iLAI
size(S⁺)


Emin,Eplu,netLW =  CanopyRTMod.RTM_diffuseS(τ_dd, ρ_dd,S⁻, S⁺, [0.0], [300.0], [0.06]);

#iLAI = canopy.LAI/canopy.nlayers;
sumLAI = [0:iLAI:nl*iLAI;];

plot(Emin',-sumLAI, label="Downwelling Thermal")
plot!(Eplu',-sumLAI, label="Upwelling Thermal")
plot!(1 .+Emin2',-sumLAI2, label="Downwelling Thermal fine")
plot!(1 .+Eplu2',-sumLAI2, label="Upwelling Thermal fine")
xlabel!("W/m2")
ylabel!("-Cum LAI")
#plot!(Eplu',1:1:21, label="Upwnwelling Thermal")

(0.98*consts.σ*290^4)


plot(netLW'-S⁺'-S⁻',-sumLAI[1:nl], label="net Thermal")
#plot!(Eplu',-sumLAI, label="Upwelling Thermal")
xlabel!("W/m2")
ylabel!("-Cum LAI")

# Hack Thermal optical properties per layer here:
nl = 200
nwl = 1
iLAI    = 0.02;#canopy.LAI/nl;
sigf = canOpt.ddf*leaf.ρ_LW + canOpt.ddb*leaf.τ_LW
sigb = canOpt.ddb*leaf.ρ_LW + canOpt.ddf*leaf.τ_LW
τ_dd = (1 .- (1-sigf)*iLAI)*ones(nwl,nl)
ρ_dd = (sigb*iLAI)*ones(nwl,nl)

1.0-τ_dd[1,1]-ρ_dd[1,1]
@show 1-exp(-sigf*iLAI)
@show (1-sigf)*iLAI

S⁺ = zeros(nwl,nl).+150*iLAI*0.98
S⁻ = zeros(nwl,nl).+150*iLAI*0.98
#S⁺[1,80:109].=200*iLAI
#S⁻[1,80:109].=200*iLAI
size(S⁺)
Emin2,Eplu2,netLW2 =  CanopyRTMod.RTM_diffuseS(τ_dd, ρ_dd,S⁻, S⁺, [0.0], [300.0], [0.06]);


sumLAI2 = [0:iLAI:nl*iLAI;];
plot(Emin2',-sumLAI2, label="Downwelling Thermal")
plot!(Eplu2',-sumLAI2, label="Upwelling Thermal")
xlabel!("W/m2")
ylabel!("-Cum LAI")

# Hack Thermal optical properties per layer here:
nl = 20
nwl = 1
iLAI    = 4/nl;0.002;#canopy.LAI/nl;
sigf = canOpt.ddf*leaf.ρ_LW + canOpt.ddb*leaf.τ_LW
sigb = canOpt.ddb*leaf.ρ_LW + canOpt.ddf*leaf.τ_LW
τ_dd = (1 - (exp(-sigf*iLAI))*ones(nwl,nl)
ρ_dd = (sigb*iLAI)*ones(nwl,nl)

#
#
#
#
# a parsing error here
# 1.0-τ_dd[1,1]-ρ_dd[1,1]
#
#
#
#
S⁺ = zeros(nwl,nl).+150*iLAI*0.98
S⁻ = zeros(nwl,nl).+150*iLAI*0.98
#S⁺[1,80:109].=200*iLAI
#S⁻[1,80:109].=200*iLAI
size(S⁺)
Emin3,Eplu3,netLW3 =  CanopyRTMod.RTM_diffuseS(τ_dd, ρ_dd,S⁻, S⁺, [0.0], [300.0], [0.06]);
sumLAI3 = [0:iLAI:nl*iLAI;];
plot(Emin2',-sumLAI2, label="Downwelling Thermal")
plot!(Eplu2',-sumLAI2, label="Upwelling Thermal")
xlabel!("W/m2")
ylabel!("-Cum LAI")

plot(Emin3',-sumLAI3, label="Downwelling Thermal")
plot!(Eplu3',-sumLAI3, label="Upwelling Thermal")
plot!(1 .+Emin2',-sumLAI2, label="Downwelling Thermal fine")
plot!(1 .+Eplu2',-sumLAI2, label="Upwelling Thermal fine")
xlabel!("W/m2")
ylabel!("-Cum LAI")
=#

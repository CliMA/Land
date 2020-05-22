# Add usual tools we use:
#using Revise
using BenchmarkTools
using Plots
pyplot()

# load Photosynthesis module:
using Land.Photosynthesis

# Specify Field Type
const FT = Float32

# Create a standard leaf with defualt parameters
leaf = leaf_params{FT}();

# Create a standard meteo structure:
met = meteo{FT}();

# Setting some standard values (dynamic-state=false forces A-Cc iterations)
leaf.dynamic_state = false


met.stab_type_stable = 2;
met.e_air = 1500;
met.T_air = 298;

# use this as the boundary layer resistance (1/gb)
ra = 0.5
modC3 = Photosynthesis.PhotoMods(
    fluorescence    = FlexasTolBerryFluorescence{FT}(),
    photosynthesis  = C3FvCBPhoto(),
    respiration     = RespirationCLM{FT}(),
    stomatal        = BallBerryStomata{FT}(g1=8),
    Jmax            = JmaxCLM{FT}(),
    Vmax            = VcmaxCLM{FT}(),
    MichaelisMenten = MM_CLM{FT}(),
    BoundaryLayer   = FixedBoundaryResistance{FT}(ra=ra),
    colimitation = CurvedColimit{FT}(0.99));

modC4 = Photosynthesis.PhotoMods(
    fluorescence    = FlexasTolBerryFluorescence{FT}(),
    photosynthesis  = C4CollatzPhoto(),
    respiration     = RespirationCLM{FT}(),
    stomatal        = BallBerryStomata{FT}(g1=4),
    Jmax            = JmaxCLM{FT}(),
    Vmax            = VcmaxCLM{FT}(),
    MichaelisMenten = MM_CLM{FT}(),
    BoundaryLayer   = FixedBoundaryResistance{FT}(ra=ra),
    colimitation = CurvedColimit{FT}(0.99));

##Again, this looks tedious and repetitive but is the easiest way to do this right now:

# Variable I want to save:
Ag_C3 = Float32[]; Ag_C4 = Float32[]
An_C3 = Float32[]; An_C4 = Float32[]
Aj_C3 = Float32[]; Aj_C4 = Float32[]
Ap_C3 = Float32[]; Ap_C4 = Float32[]
Ac_C3 = Float32[]; Ac_C4 = Float32[]
Cc_C3 = Float32[]; Cc_C4 = Float32[]
gs_C3 = Float32[]; gs_C4 = Float32[]

# Set Ca to 400ppm
met.Ca = 400

# Set leaf temperature:
leaf.T = 298.0

# Specify Vcmax
leaf.Vcmax25 = 90
leaf.Jmax25 = 90*1.9

# Specify curvature for J
leaf.θ_j = 0.995

APAR = collect(0:10:1700)
for leaf.APAR in APAR
    # Run C3 photosynthesis:
    LeafPhotosynthesis!(modC3, leaf, met);
    # Save leaf variables:
    push!(An_C3, leaf.An);push!(Ag_C3, leaf.Ag);push!(Aj_C3, leaf.Aj);push!(Ap_C3, leaf.Ap); push!(Ac_C3, leaf.Ac); push!(Cc_C3, leaf.Cc); push!(gs_C3, leaf.gs)

    # Run C4 photosynthesis:
    LeafPhotosynthesis!(modC4, leaf, met);
    # Save leaf variables:
    push!(An_C4, leaf.An);push!(Ag_C4, leaf.Ag);push!(Aj_C4, leaf.Aj);push!(Ap_C4, leaf.Ap); push!(Ac_C4, leaf.Ac); push!(Cc_C4, leaf.Cc);; push!(gs_C4, leaf.gs)
end

# Testing some times, how long does this take (as we need to run it globally, it has to be efficient)?
@btime LeafPhotosynthesis!(modC3, leaf, met);
@btime LeafPhotosynthesis!(modC4, leaf, met);

##plot(APAR, An,  label="An")
plot(APAR, Ag_C3,color=[:black],lw=1, alpha=0.7, label="Ag C3")
plot!(APAR, Ac_C3,line=(:dash,2), label="Ac C3")
plot!(APAR, Aj_C3,line=(:dash,2), label="Aj C3")
plot!(APAR, Ap_C3,line=(:dash,2), label="Ap C3" )
xlabel!("APAR [μmol/m2/s]")
ylabel!("Aᵢ [μmol/m2/s]")
title!("C3 photosynthesis light response")
##plot!(APAR, An2)

##plot(APAR, An,  label="An")
plot(APAR, Ag_C4,color=[:black],lw=1, alpha=0.7, label="Ag C4")
plot!(APAR, Ac_C4,line=(:dash,2), label="Ac C4")
plot!(APAR, Aj_C4,line=(:dash,2), label="Aj C4")
plot!(APAR, Ap_C4,line=(:dash,2), label="Ap C4" )
xlabel!("APAR [μmol/m2/s]")
ylabel!("Aᵢ [μmol/m2/s]")
title!("C4 photosynthesis light response")

plot(APAR, Cc_C4/met.Ca,color=[:black],lw=2, alpha=0.7, label="Cc/Ca C4")
plot!(APAR, Cc_C3/met.Ca,color=[:orange],lw=2, alpha=0.7, label="Cc/Ca C3")
xlabel!("APAR [μmol/m²/s]")
ylabel!("Cc/Ca [-]")

plot(APAR, gs_C4,color=[:black],lw=2, alpha=0.7, label="gs C4")
plot!(APAR, gs_C3,color=[:orange],lw=2, alpha=0.7, label="gs C3")
xlabel!("APAR [μmol/m²/s]")
ylabel!("gs")

##Again, this looks tedious and repetitive but is the easiest way to do this right now:

# Now I want to vary T, APAR and CO2:
APAR = [100.0, 250.0, 500.0, 1000.0, 1500.0]
CO2  = collect(2:2:800)
T    = collect(265:1:310)

n1 = length(APAR);
n2 = length(CO2);
n3 = length(T);

Ag_C3 = zeros(n1,n2,n3); Ag_C4 = zeros(n1,n2,n3)
An_C3 = zeros(n1,n2,n3); An_C4 = zeros(n1,n2,n3)
Aj_C3 = zeros(n1,n2,n3); Aj_C4 = zeros(n1,n2,n3)
Ap_C3 = zeros(n1,n2,n3); Ap_C4 = zeros(n1,n2,n3)
Ac_C3 = zeros(n1,n2,n3); Ac_C4 = zeros(n1,n2,n3)
Cc_C3 = zeros(n1,n2,n3); Cc_C4 = zeros(n1,n2,n3)
gs_C3 = zeros(n1,n2,n3); gs_C4 = zeros(n1,n2,n3)

# Set Ca to 400ppm
met.Ca = 400

# Set leaf temperature:
leaf.T = 298.0

# Specify Vcmax
leaf.Vcmax25 = 90
leaf.Jmax25 = 90*1.9

# Specify curvature for J
leaf.θ_j = 0.995

# Run this over all potential 3D dimensions:

# I really like the compact form of nested loops in Julia!
for iA in eachindex(APAR), iC in eachindex(CO2), iT in eachindex(T)
    met.Ca    = CO2[iC];
    # make sure we have a reasonable start:
    leaf.Cc   = met.Ca;
    leaf.Cs   = met.Ca;
    leaf.gs   = 0.1;
    leaf.T    = T[iT];
    leaf.APAR = APAR[iA];

    # Run C3 photosynthesis:
    LeafPhotosynthesis!(modC3, leaf, met);

    # Save leaf variables:
    An_C3[iA,iC,iT]=leaf.An;
    Ag_C3[iA,iC,iT]=leaf.Ag;
    Aj_C3[iA,iC,iT]=leaf.Aj;
    Ap_C3[iA,iC,iT]=leaf.Ap;
    Ac_C3[iA,iC,iT]=leaf.Ac;
    Cc_C3[iA,iC,iT]=leaf.Cc;
    gs_C3[iA,iC,iT]=leaf.gs;

    # Run C4 photosynthesis:
    leaf.Cc   = met.Ca;
    leaf.Cs   = met.Ca;
    leaf.gs   = 0.1;
    LeafPhotosynthesis!(modC4, leaf, met);

    # Save leaf variables:
    An_C4[iA,iC,iT]=leaf.An;
    Ag_C4[iA,iC,iT]=leaf.Ag;
    Aj_C4[iA,iC,iT]=leaf.Aj;
    Ap_C4[iA,iC,iT]=leaf.Ap;
    Ac_C4[iA,iC,iT]=leaf.Ac;
    Cc_C4[iA,iC,iT]=leaf.Cc;
    gs_C4[iA,iC,iT]=leaf.gs;
end

##Let's take one slice in APAR space:
i = 4

# and plot:
contourf(T.-273.15, CO2, An_C3[i,:,:])
xlabel!("T [°C]")
ylabel!("Ambient CO₂ [ppm]")
title!("C3 , An [μmol/m²/s] at APAR=$(APAR[i])")

# Same for C4 plants, why is it so different??

contourf(T.-273.15, CO2, An_C4[i,:,:])
xlabel!("T [°C]")
ylabel!("Ambient CO₂ [ppm]")
title!("C4 , An [μmol/m²/s] at APAR=$(APAR[i])")

iA = 4; iT=31
plot( CO2, Ag_C3[iA,:,iT],color=[:black],lw=2, alpha=0.7, label="Ag C3")
plot!(CO2, Ac_C3[iA,:,iT],line=(:dash,2), label="Ac C3")
plot!(CO2, Aj_C3[iA,:,iT],line=(:dash,2), label="Aj C3")
plot!(CO2, Ap_C3[iA,:,iT],line=(:dash,2), label="Ap C3" )
xlabel!("CO₂ [ppm]")
ylabel!("Aᵢ [μmol/m2/s]")
title!("Ambient C3 CO₂ response, T=$(T[iT]-273), APAR=$(APAR[iA])")

iA = 4; iT=31
plot(CO2, Ag_C4[iA,:,iT],color=[:black],lw=2, alpha=0.7, label="Ag C4")
plot!(CO2, Ac_C4[iA,:,iT],line=(:dash,2), label="Ac C4")
plot!(CO2, Aj_C4[iA,:,iT],line=(:dash,2), label="Aj C4")
plot!(CO2, Ap_C4[iA,:,iT],line=(:dash,2), label="Ap C4" )
xlabel!("CO₂ [ppm]")
ylabel!("Aᵢ [μmol/m2/s]")
title!("Ambient C4 CO₂ response, T=$(T[iT]-273), APAR=$(APAR[iA])")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl


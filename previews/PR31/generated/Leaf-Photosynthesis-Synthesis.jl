# Add usual tools we use:
##using Revise
using BenchmarkTools
using PyPlot

# load Photosynthesis module:
using Land.Photosynthesis

# Specify Field Type
const FT = Float64

# Create a standard leaf with defualt parameters
leaf3 = Leaf{FT}(APAR=1200, Vcmax25=90, Jmax25=90*1.9, Vpmax25=100, Rd25=1);
leaf4 = Leaf{FT}(APAR=1200, Vcmax25=90, Jmax25=90*1.9, Vpmax25=100, Rd25=1);
# Create a standard meteo structure:
envir  = AirLayer{FT}();

# Setting some standard values (dynamic-state=false forces A-Cc iterations)
##leaf.dynamic_state = false

# use this as the boundary layer resistance (1/gb)
# deprecated as there are g_bw and g_bc terms in Leaf struct
# ra = FT(0.5)

# C3 Photosynthesis

modC3 = C3CLM(FT)
modC3.Col = Photosynthesis.CurvedColimit{FT}();
modC3.Sto = Photosynthesis.ESMBallBerry{FT}(g1 = 16)

# C4 Photosynthesis
modC4 = C4CLM(FT)
modC4.Col = Photosynthesis.CurvedColimit{FT}();
modC4.Sto = Photosynthesis.ESMBallBerry{FT}(g1 = 8)

##Again, this looks tedious and repetitive but is the easiest way to do this right now:

# Variable I want to save:
Ag_C3 = Float32[]; Ag_C4 = Float32[];
An_C3 = Float32[]; An_C4 = Float32[];
Aj_C3 = Float32[]; Aj_C4 = Float32[];
Ap_C3 = Float32[]; Ap_C4 = Float32[];
Ac_C3 = Float32[]; Ac_C4 = Float32[];
Cc_C3 = Float32[]; Cc_C4 = Float32[];
gs_C3 = Float32[]; gs_C4 = Float32[];

APAR = collect(FT, 0:10:1700)
for _APAR in APAR
    leaf3.APAR = _APAR;
    leaf4.APAR = _APAR;
    # Run C3 photosynthesis (general photosynthesis model):
    leaf_photo_from_envir!(modC3, leaf3, envir, modC3.Sto);
    # Save leaf variables:
    push!(An_C3, leaf3.An); push!(Ag_C3, leaf3.Ag);
    push!(Aj_C3, leaf3.Aj); push!(Ap_C3, leaf3.Ap);
    push!(Ac_C3, leaf3.Ac); push!(Cc_C3, leaf3.p_i);
    push!(gs_C3, leaf3.g_sw);

    # Run C4 photosynthesis:
    leaf_photo_from_envir!(modC4, leaf4, envir, modC4.Sto);
    # Save leaf variables:
    push!(An_C4, leaf4.An); push!(Ag_C4, leaf4.Ag);
    push!(Aj_C4, leaf4.Aj); push!(Ap_C4, leaf4.Ap);
    push!(Ac_C4, leaf4.Ac); push!(Cc_C4, leaf4.p_i);
    push!(gs_C4, leaf4.g_sw);
end

# Testing some times, how long does this take (as we need to run it globally, it has to be efficient)?
# Slow for now, because of unnecessary allocations
# Will improve when the structs are cleaned up
#@btime leaf_photo_from_envir!(modC3, leaf3, envir, modC3.Sto);
#@btime leaf_photo_from_envir!(modC4, leaf4, envir, modC4.Sto);

##plot(APAR, An,  label="An")
figure()
plot(APAR, Ag_C3,color=:black,lw=2, alpha=0.7, label="Ag C3")
plot(APAR, Ac_C3, ls="--", lw=2, label="Ac C3")
plot(APAR, Aj_C3, ls="--", lw=2, label="Aj C3")
plot(APAR, Ap_C3, ls="--", lw=2, label="Ap C3" )
xlabel("APAR [μmol/m2/s]")
ylabel("Aᵢ [μmol/m2/s]")
title("C3 photosynthesis light response")
legend()
gcf()

##plot(APAR, An,  label="An")
figure()
plot(APAR, Ag_C4,color=:black, lw=2, alpha=0.7, label="Ag C4")
plot(APAR, Ac_C4, ls="--", lw=2, label="Ac C4")
plot(APAR, Aj_C4, ls="--", lw=2, label="Aj C4")
plot(APAR, Ap_C4, ls="--", lw=2, label="Ap C4" )
xlabel("APAR [μmol/m2/s]")
ylabel("Aᵢ [μmol/m2/s]")
title("C4 photosynthesis light response")
legend()
gcf()

figure()
plot(APAR, Cc_C4/envir.p_a,color=:black ,lw=2, alpha=0.7, label="Cc/Ca C4")
plot(APAR, Cc_C3/envir.p_a,color=:orange,lw=2, alpha=0.7, label="Cc/Ca C3")
xlabel("APAR [μmol/m²/s]")
ylabel("Cc/Ca [-]")
legend()
gcf()

#This part has been broken by Yuije, please fix it
figure()
plot(APAR, gs_C4,color=:black ,lw=2, alpha=0.7, label="gs C4")
plot(APAR, gs_C3,color=:orange,lw=2, alpha=0.7, label="gs C3")
xlabel("APAR [μmol/m²/s]")
ylabel("gs")
legend()
gcf()

##Again, this looks tedious and repetitive but is the easiest way to do this right now:

# Now I want to vary T, APAR and CO2:
# Start CO2 from 100 ppm to make sure it is higher than Γ*
APAR = [100.0, 250.0, 500.0, 1000.0, 1500.0]
CO2  = collect(100:20:800)
T    = collect(274:2:310)

n1 = length(APAR);
n2 = length(CO2);
n3 = length(T);

Ag_C3 = zeros(n1,n2,n3); Ag_C4 = zeros(n1,n2,n3);
An_C3 = zeros(n1,n2,n3); An_C4 = zeros(n1,n2,n3);
Aj_C3 = zeros(n1,n2,n3); Aj_C4 = zeros(n1,n2,n3);
Ap_C3 = zeros(n1,n2,n3); Ap_C4 = zeros(n1,n2,n3);
Ac_C3 = zeros(n1,n2,n3); Ac_C4 = zeros(n1,n2,n3);
Cc_C3 = zeros(n1,n2,n3); Cc_C4 = zeros(n1,n2,n3);
gs_C3 = zeros(n1,n2,n3); gs_C4 = zeros(n1,n2,n3);

# Run this over all potential 3D dimensions:

# I really like the compact form of nested loops in Julia!
for iA in eachindex(APAR), iC in eachindex(CO2), iT in eachindex(T)
    #println(iA, "/", iC, "/", iT)
    envir.p_a  = CO2[iC]/10;
    leaf3.T    = T[iT];
    leaf3.APAR = APAR[iA];
    leaf4.T    = T[iT];
    leaf4.APAR = APAR[iA];

    # Run C3 photosynthesis:
    leaf_photo_from_envir!(modC3, leaf3, envir, modC3.Sto);

    # Save leaf variables:
    An_C3[iA,iC,iT]=leaf3.An;
    Ag_C3[iA,iC,iT]=leaf3.Ag;
    Aj_C3[iA,iC,iT]=leaf3.Aj;
    Ap_C3[iA,iC,iT]=leaf3.Ap;
    Ac_C3[iA,iC,iT]=leaf3.Ac;
    Cc_C3[iA,iC,iT]=leaf3.p_i;
    gs_C3[iA,iC,iT]=leaf3.g_sw;

    # Run C4 photosynthesis:
    leaf_photo_from_envir!(modC4, leaf4, envir, modC4.Sto);

    # Save leaf variables:
    An_C4[iA,iC,iT]=leaf4.An;
    Ag_C4[iA,iC,iT]=leaf4.Ag;
    Aj_C4[iA,iC,iT]=leaf4.Aj;
    Ap_C4[iA,iC,iT]=leaf4.Ap;
    Ac_C4[iA,iC,iT]=leaf4.Ac;
    Cc_C4[iA,iC,iT]=leaf4.p_i;
    gs_C4[iA,iC,iT]=leaf4.g_sw;
end

##Let's take one slice in APAR space:
i = 4

# and plot:
figure()
contourf(T.-273.15, CO2, An_C3[i,:,:])
xlabel("T [°C]")
ylabel("Ambient CO₂ [ppm]")
title("C3 , An [μmol/m²/s] at APAR=$(APAR[i])")
colorbar()
gcf()

# Same for C4 plants, why is it so different??

figure()
contourf(T.-273.15, CO2, Cc_C4[i,:,:])
xlabel("T [°C]")
ylabel("Ambient CO₂ [ppm]")
title("C4 , Cc [Pa] at APAR=$(APAR[i])")
colorbar()
gcf()

# Same for C4 plants, why is it so different??

figure()
contourf(T.-273.15, CO2, An_C4[i,:,:])
xlabel("T [°C]")
ylabel("Ambient CO₂ [ppm]")
title("C4 , An [μmol/m²/s] at APAR=$(APAR[i])")
colorbar()
gcf()

iA = 4; iT=12
figure()
plot( CO2, Ag_C3[iA,:,iT],color=:black,lw=2, alpha=0.7, label="Ag C3")
plot(CO2, Ac_C3[iA,:,iT], ls="--", lw=2, label="Ac C3")
plot(CO2, Aj_C3[iA,:,iT], ls="--", lw=2, label="Aj C3")
plot(CO2, Ap_C3[iA,:,iT], ls="--", lw=2, label="Ap C3" )
xlabel("CO₂ [ppm]")
ylabel("Aᵢ [μmol/m2/s]")
title("Ambient C3 CO₂ response, T=$(T[iT]-273), APAR=$(APAR[iA])")
legend()
gcf()

iA = 4; iT=12
figure()
plot(CO2, Ag_C4[iA,:,iT],color=:black,lw=2, alpha=0.7, label="Ag C4")
plot(CO2, Ac_C4[iA,:,iT], ls="--", lw=2, label="Ac C4")
plot(CO2, Aj_C4[iA,:,iT], ls="--", lw=2, label="Aj C4")
plot(CO2, Ap_C4[iA,:,iT], ls="--", lw=2, label="Ap C4")
xlabel("CO₂ [ppm]")
ylabel("Aᵢ [μmol/m2/s]")
title("Ambient C4 CO₂ response, T=$(T[iT]-273), APAR=$(APAR[iA])")
legend()
gcf()

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl


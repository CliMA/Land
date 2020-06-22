# Add usual tools we use:
#using Revise
using BenchmarkTools

using PyPlot

# load Photosynthesis module:
using Land.Photosynthesis

# Specify Field Type
const FT = Float32

# This looks a bit more tedious here than it needs to but in reality

vcmax_Bernacchi = Float32[]
vcmax_CLM = Float32[]

# Define T-range
Tleaf = collect(FT,260.0:1.0:315.0)

# Run through temperatures and save Vcmax values:
td_vc_bernacchi = Photosynthesis.VcmaxTDBernacchi(FT)
td_vc_clm       = Photosynthesis.VcmaxTDCLM(FT)
for T in Tleaf
    _Vcmax = Photosynthesis.photo_TD_from_val(td_vc_bernacchi, FT(100.0), T)
    push!(vcmax_Bernacchi, _Vcmax)
    _Vcmax = Photosynthesis.photo_TD_from_val(td_vc_clm, FT(100.0), T)
    push!(vcmax_CLM, _Vcmax)
end

figure()
plot(Tleaf .- 273.15, vcmax_Bernacchi, label="Bernacchi", lw=2)
plot(Tleaf .- 273.15, vcmax_CLM,  label="CLM5", lw=2)
ylabel("Vcmax (µmol/m²/s)")
xlabel("Leaf Temperature (°C)")
legend()
gcf()

# Here, we only have one implementation:
Γ_CLM = Float32[]
Tleaf = collect(FT,260.0:1.0:315.0)
td_gamma_clm = Photosynthesis.ΓStarTDCLM(FT)

for T in Tleaf
    _ΓStar = Photosynthesis.photo_TD_from_set(td_gamma_clm, T)
    push!(Γ_CLM, _ΓStar)
end

figure()
plot(Tleaf .- 273.15, Γ_CLM, label="Γ_CLM")
ylabel("Γ⋆ (Pa)")
xlabel("Leaf Temperature (°C)")
legend()
gcf()

Jmax_Bernacchi = Float32[]
Jmax_CLM = Float32[]
Tleaf = collect(FT,260.0:1.0:315.0)

td_j_bernacchi = Photosynthesis.JmaxTDBernacchi(FT)
td_j_clm = Photosynthesis.JmaxTDCLM(FT)

for T in Tleaf
    _Jmax = Photosynthesis.photo_TD_from_val(td_j_bernacchi, FT(100.0), T)
    push!(Jmax_Bernacchi, _Jmax)
    _Jmax = Photosynthesis.photo_TD_from_val(td_j_clm, FT(100.0), T)
    push!(Jmax_CLM, _Jmax)
end

figure()
plot(Tleaf .- 273.15, Jmax_Bernacchi, label="Bernacchi")
plot(Tleaf .- 273.15, Jmax_CLM,  label="CLM5")
ylabel("Jmax (µmol/m2/s)")
xlabel("Leaf Temperature (°C)")
legend()
gcf()

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl


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
for T in Tleaf
    _Vcmax = Photosynthesis.get_vmax( Photosynthesis.VcmaxTDBernacchi(FT), FT(100.0), T )
    push!(vcmax_Bernacchi, _Vcmax)
    _Vcmax = Photosynthesis.get_vmax( Photosynthesis.VcmaxTDCLM(FT), FT(100.0), T )
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

for T in Tleaf
    _ΓStar = Photosynthesis.get_Γ_star( Photosynthesis.ΓStarTDCLM(FT), T )
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

for T in Tleaf
    _Jmax = Photosynthesis.get_jmax( Photosynthesis.JmaxTDBernacchi(FT), FT(100.0), T )
    push!(Jmax_Bernacchi, _Jmax)
    _Jmax = Photosynthesis.get_jmax( Photosynthesis.JmaxTDCLM(FT), FT(100.0), T )
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


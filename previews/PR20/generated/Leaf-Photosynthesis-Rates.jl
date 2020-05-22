# Add usual tools we use:
using Revise
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

fieldnames(Photosynthesis.leaf_params)

# This looks a bit more tedious here than it needs to but in reality
vcmax_Bernacchi = Float32[];vcmax_CLM = Float32[]

# Define T-range
Tleaf = 260:1:315

# Run through temperatures and save Vcmax values:
for T in Tleaf
    leaf.T = T
    Photosynthesis.max_carboxylation_rate!(VcmaxBernacchi{Float32}(), leaf)
    push!(vcmax_Bernacchi, leaf.Vcmax)
    Photosynthesis.max_carboxylation_rate!(VcmaxCLM{Float32}(), leaf)
    push!(vcmax_CLM, leaf.Vcmax)
end

plot(Tleaf .- 273.15, vcmax_Bernacchi, label="Bernacchi", lw=2)
plot!(Tleaf .- 273.15, vcmax_CLM,  label="CLM5", lw=2)
ylabel!("Vcmax (µmol/m²/s)")
xlabel!("Leaf Temperature (°C)")

# Here, we only have one implementation:
Γ_CLM = Float32[]
Tleaf = 260:1:305

for T in Tleaf
    leaf.T = T
    Photosynthesis.michaelis_menten_constants!(MM_CLM(), leaf)
    push!(Γ_CLM, leaf.Γstar)
end

plot(Tleaf .- 273.15, Γ_CLM, label="Γ_CLM")
ylabel!("Γ⋆ (Pa)")
xlabel!("Leaf Temperature (C)")

Jmax_Bernacchi = Float32[]
Jmax_CLM = Float32[]

Tleaf = 260:1:315

for T in Tleaf
    leaf.T = T
    Photosynthesis.max_electron_transport_rate!(JmaxBernacchi{Float32}(), leaf)
    push!(Jmax_Bernacchi, leaf.Jmax)
    Photosynthesis.max_electron_transport_rate!(JmaxCLM{Float32}(), leaf)
    push!(Jmax_CLM, leaf.Jmax)
end

plot(Tleaf .- 273.15, Jmax_Bernacchi, label="Bernacchi")
plot!(Tleaf .- 273.15, Jmax_CLM,  label="CLM5")
ylabel!("Jmax (µmol/m2/s)")
xlabel!("Leaf Temperature (C)")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl



# # Rate Constants, T-dependence
# Here, we will just briefly summarize how temperature affects enzymatic rate constants and how much uncertainty there actually is in literature (in terms of how to best define them and how plastic some of the variables are)

## Add usual tools we use:
using Revise
using BenchmarkTools
using Plots
pyplot()
#----------------------------------------------------------------------------

## load Photosynthesis module: 
using Land.Photosynthesis
#----------------------------------------------------------------------------

## Specify Field Type
const FT = Float32
#----------------------------------------------------------------------------

## Create a standard leaf with defualt parameters
leaf = leaf_params{FT}();

## Create a standard meteo structure:
met = meteo{FT}();
#----------------------------------------------------------------------------

fieldnames(Photosynthesis.leaf_params)
#----------------------------------------------------------------------------

# ## T-dependence of V$_{c,max}$
# In literature, there are different implementations of the temperature dependence of enzymatic reaction rates. Common among all of them is a typical Arrhenius formulation with Activation energy $E_a$, which leads to a temperature dependence of a quantity $V$ given the standard value defined at 25°C as $T_{ref}$.
# 
# $$V(T)= V(T_{ref})\, \underbrace{\exp\left(\frac{E_a}{RT_{ref}}-\frac{E_a}{RT}\right)}_{\text{Activation}}$$
# 
# Other formulations add de-activation of proteins due to denaturalization at higher temperatures:
# $$V(T)= V(T_{ref})\, \exp\left(\frac{E_a}{RT_{ref}}-\frac{E_a}{RT}\right) \underbrace{\frac{1+\exp\left((\Delta ST_{ref}-H_d)/RT_{ref}\right)}{1+\exp\left((\Delta ST-H_d)/RT\right)}}_{\text{de-activation}} $$
# which includes an entropy term $\Delta S$ and the energy for deactivation $H_d$.
# 
# To illustrate the differences, we show two different implementations, one only using the activation part based on Bernacchi et al 2001[^1] and one with the typical CLM5 implementation.
# 
# [^1]: Bernacchi, C.J., Singsaas, E.L., Pimentel, C., Portis Jr, A.R. and Long, S.P., 2001. Improved temperature response functions for models of Rubisco‐limited photosynthesis. Plant, Cell & Environment, 24(2), pp.253-259.
# 

## This looks a bit more tedious here than it needs to but in reality
vcmax_Bernacchi = Float32[];vcmax_CLM = Float32[]

## Define T-range
Tleaf = 260:1:315

## Run through temperatures and save Vcmax values:
for T in Tleaf
    leaf.T = T
    Photosynthesis.max_carboxylation_rate!(VcmaxBernacchi{Float32}(), leaf)
    push!(vcmax_Bernacchi, leaf.Vcmax)
    Photosynthesis.max_carboxylation_rate!(VcmaxCLM{Float32}(), leaf)
    push!(vcmax_CLM, leaf.Vcmax)
end
#----------------------------------------------------------------------------

plot(Tleaf .- 273.15, vcmax_Bernacchi, label="Bernacchi", lw=2)
plot!(Tleaf .- 273.15, vcmax_CLM,  label="CLM5", lw=2)
ylabel!("Vcmax (µmol/m²/s)")
xlabel!("Leaf Temperature (°C)")
#----------------------------------------------------------------------------

# ### CO$_2$ compensation point $\Gamma_\star$

## Here, we only have one implementation:
Γ_CLM = Float32[]
Tleaf = 260:1:305

for T in Tleaf
    leaf.T = T
    Photosynthesis.michaelis_menten_constants!(MM_CLM(), leaf)
    push!(Γ_CLM, leaf.Γstar)
end
#----------------------------------------------------------------------------

plot(Tleaf .- 273.15, Γ_CLM, label="Γ_CLM")
ylabel!("Γ⋆ (Pa)")
xlabel!("Leaf Temperature (C)")
#----------------------------------------------------------------------------

# ## T-dependence of J$_{max}$

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
#----------------------------------------------------------------------------

plot(Tleaf .- 273.15, Jmax_Bernacchi, label="Bernacchi")
plot!(Tleaf .- 273.15, Jmax_CLM,  label="CLM5")
ylabel!("Jmax (µmol/m2/s)")
xlabel!("Leaf Temperature (C)")
#----------------------------------------------------------------------------

# --- 
# ## Summary 
# The way land surface models implement temperature variations in V$_{c,max}$ is still rather variable and there is no clear consensus in what model is best. In addition, gas exchange datasets performed at different temperatures might actually alias different confounding factors into the V$_{c,max}$ determination, so that the actually fitted V$_{c,max}$ strongly depends on the model formulation. E.g. mesophyll conductance is usually ignore in deriving V$_{c,max}$ but might be temperature dependent, which can cause T-dependent errors in V$_{c,max}$.
# 
# That said, having the option to use different model forumations and actually using the ones that are implemented in land surface models directly when fitting leaf level measurements seems to be a prudent way to improve models in the future. More leaf level data might help in the future, especially looking at temperature dependencies, mesophyll conductance, and fluorescence yields. So many chemical and phsyical aspects change with temperature that it is hard to isolate impacting factors. Now, at least, we can run the models directly in the REPL, which is great for prototyping anc validating the model on the individual module level without any effort in terms of additional coding. 

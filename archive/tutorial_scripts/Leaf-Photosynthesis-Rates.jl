#=
# # Rate Constants, T-dependence
# Here, we will just briefly summarize how temperature affects enzymatic rate constants and how much uncertainty there actually is in literature (in terms of how to best define them and how plastic some of the variables are)

## Add usual tools we use:
#using Revise
using BenchmarkTools

# Use PyPlot to plot figures
using PyPlot

#----------------------------------------------------------------------------

## load Photosynthesis module:
using Land.Photosynthesis
#----------------------------------------------------------------------------

## Specify Field Type
const FT = Float32
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

vcmax_Bernacchi = Float32[]
vcmax_CLM = Float32[]

## Define T-range
Tleaf = collect(FT,260.0:1.0:315.0)

## Run through temperatures and save Vcmax values:
td_vc_bernacchi = Photosynthesis.VcmaxTDBernacchi(FT)
td_vc_clm       = Photosynthesis.VcmaxTDCLM(FT)
for T in Tleaf
    _Vcmax = Photosynthesis.photo_TD_from_val(td_vc_bernacchi, FT(100.0), T)
    push!(vcmax_Bernacchi, _Vcmax)
    _Vcmax = Photosynthesis.photo_TD_from_val(td_vc_clm, FT(100.0), T)
    push!(vcmax_CLM, _Vcmax)
end
#----------------------------------------------------------------------------

figure()
plot(Tleaf .- 273.15, vcmax_Bernacchi, label="Bernacchi", lw=2)
plot(Tleaf .- 273.15, vcmax_CLM,  label="CLM5", lw=2)
ylabel("Vcmax (µmol/m²/s)")
xlabel("Leaf Temperature (°C)")
legend()
gcf()
#----------------------------------------------------------------------------

# ### CO$_2$ compensation point $\Gamma_\star$

## Here, we only have one implementation:
Γ_CLM = Float32[]
Tleaf = collect(FT,260.0:1.0:315.0)
td_gamma_clm = Photosynthesis.ΓStarTDCLM(FT)

for T in Tleaf
    _ΓStar = Photosynthesis.photo_TD_from_set(td_gamma_clm, T)
    push!(Γ_CLM, _ΓStar)
end
#----------------------------------------------------------------------------
figure()
plot(Tleaf .- 273.15, Γ_CLM, label="Γ_CLM")
ylabel("Γ⋆ (Pa)")
xlabel("Leaf Temperature (°C)")
legend()
gcf()
#----------------------------------------------------------------------------

# ## T-dependence of J$_{max}$

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
#----------------------------------------------------------------------------
figure()
plot(Tleaf .- 273.15, Jmax_Bernacchi, label="Bernacchi")
plot(Tleaf .- 273.15, Jmax_CLM,  label="CLM5")
ylabel("Jmax (µmol/m2/s)")
xlabel("Leaf Temperature (°C)")
legend()
gcf()
#----------------------------------------------------------------------------

# ---
# ## Summary
# The way land surface models implement temperature variations in V$_{c,max}$ is still rather variable and there is no clear consensus in what model is best. In addition, gas exchange datasets performed at different temperatures might actually alias different confounding factors into the V$_{c,max}$ determination, so that the actually fitted V$_{c,max}$ strongly depends on the model formulation. E.g. mesophyll conductance is usually ignore in deriving V$_{c,max}$ but might be temperature dependent, which can cause T-dependent errors in V$_{c,max}$.
#
# That said, having the option to use different model forumations and actually using the ones that are implemented in land surface models directly when fitting leaf level measurements seems to be a prudent way to improve models in the future. More leaf level data might help in the future, especially looking at temperature dependencies, mesophyll conductance, and fluorescence yields. So many chemical and phsyical aspects change with temperature that it is hard to isolate impacting factors. Now, at least, we can run the models directly in the REPL, which is great for prototyping anc validating the model on the individual module level without any effort in terms of additional coding.
=#

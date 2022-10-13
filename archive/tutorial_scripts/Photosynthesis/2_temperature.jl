# # Rate Constants and T-Dependency
# 
# Here, we will just briefly summarize how temperature affects enzymatic rate constants and how much uncertainty there actually is in literature (in terms of how to best define them and how plastic some of the variables are)

## Add usual tools we use:
using BenchmarkTools
using PyPlot

#----------------------------------------------------------------------------

## load Photosynthesis module: 
using Land.Photosynthesis
#----------------------------------------------------------------------------

## Specify Field Type
const FT = Float32;
#----------------------------------------------------------------------------

# ## T-dependence
# In literature, there are different implementations of the temperature dependence of enzymatic reaction rates. Common among all of them is a typical Arrhenius formulation with Activation energy $E_a$, which leads to a temperature dependence of a quantity $V$ given the standard value defined at 25 °C as $T_{ref}$.
# 
# $$V(T)= V(T_{ref})\, \underbrace{\exp\left(\frac{E_a}{RT_{ref}}-\frac{E_a}{RT}\right)}_{\text{Activation}}$$
# 
# Other formulations add de-activation of proteins due to denaturalization at higher temperatures:
# 
# $$V(T)= V(T_{ref})\, \exp\left(\frac{E_a}{RT_{ref}}-\frac{E_a}{RT}\right) \underbrace{\frac{1+\exp\left[(\Delta ST_{ref}-H_d)/RT_{ref}\right]}{1+\exp\left[(\Delta ST-H_d)/RT\right]}}_{\text{de-activation}}$$
# 
# which includes an entropy term $\Delta S$ and the energy for deactivation $H_d$.
# 
# To illustrate the differences, we show two different implementations, one only using the activation part based on Bernacchi et al 2001[^1] and one with the typical CLM5 implementation.
# 
# The T-Dependency applies to $V_{max}$, $J_{max}$, $K_{c}$, $K_{o}$, $K_{pep}$, Respiration, and $\Gamma^{*}$, and there are multiple parameter sets for each of the parameters. Below, we show the differences of the sets for each of $V_{max}$, $J_{max}$, $K_{c}$, $K_{o}$, $K_{pep}$, Respiration, and $\Gamma^{*}$.
# 
# [^1]: Bernacchi, C.J., Singsaas, E.L., Pimentel, C., Portis Jr, A.R. and Long, S.P., 2001. Improved temperature response functions for models of Rubisco‐limited photosynthesis. Plant, Cell & Environment, 24(2), pp.253-259.
# 
# ### $V_{cmax}$
#

## Bernacchi set
vc_td_bernacchi = Photosynthesis.VcmaxTDBernacchi(FT)
## CLM set
vc_td_clm       = Photosynthesis.VcmaxTDCLM(FT)
## Leuning set
vc_td_leuning   = Photosynthesis.VcmaxTDLeuning(FT)

## an array of temperature from 10 to 50 degree C
t_array = collect(FT, 10:50) .+ FT(273.15)
vcm_bernacchi = photo_TD_from_val(vc_td_bernacchi, FT(60), t_array)
vcm_clm       = photo_TD_from_val(vc_td_clm      , FT(60), t_array)
vcm_leuning   = photo_TD_from_val(vc_td_leuning  , FT(60), t_array)

## plot the TD curves
figure(figsize=(8,6), dpi=100)
tight_layout(true)
plot(t_array, vcm_bernacchi, "r-", label="Bernacchi")
plot(t_array, vcm_clm      , "g-", label="CLM5"     )
plot(t_array, vcm_leuning  , "b-", label="Leuning"  )
xlabel("Temperature (K)"            , fontsize=16)
ylabel("\$V_{cmax}\$ (μmol m⁻² s⁻¹)", fontsize=16)
legend(loc="upper left")
gcf()
#----------------------------------------------------------------------------

# ### $V_{pmax}$
#

## Boyd set
vp_td_boyd = Photosynthesis.VpmaxTDBoyd(FT)

## an array of temperature from 10 to 50 degree C
vpm_boyd = photo_TD_from_val(vp_td_boyd, FT(60), t_array)

## plot the TD curves
figure(figsize=(8,6), dpi=100)
tight_layout(true)
plot(t_array, vpm_boyd, "r-", label="Boyd")
xlabel("Temperature (K)"            , fontsize=16)
ylabel("\$V_{pmax}\$ (μmol m⁻² s⁻¹)", fontsize=16)
legend(loc="upper left")
gcf()
#----------------------------------------------------------------------------

# ### $J_{max}$
#

## Bernacchi set
j_td_bernacchi = Photosynthesis.JmaxTDBernacchi(FT)
## CLM set
j_td_clm       = Photosynthesis.JmaxTDCLM(FT)
## Leuning set
j_td_leuning   = Photosynthesis.JmaxTDLeuning(FT)

## an array of temperature from 10 to 50 degree C
jm_bernacchi = photo_TD_from_val(j_td_bernacchi, FT(120), t_array)
jm_clm       = photo_TD_from_val(j_td_clm      , FT(120),  t_array)
jm_leuning   = photo_TD_from_val(j_td_leuning  , FT(120),  t_array)

## plot the TD curves
figure(figsize=(8,6), dpi=100)
tight_layout(true)
plot(t_array, jm_bernacchi, "r-", label="Bernacchi")
plot(t_array, jm_clm      , "g-", label="CLM5"     )
plot(t_array, jm_leuning  , "b-", label="Leuning"  )
xlabel("Temperature (K)"           , fontsize=16)
ylabel("\$J_{max}\$ (μmol m⁻² s⁻¹)", fontsize=16)
legend(loc="upper left")
gcf()
#----------------------------------------------------------------------------

# ### $K_{c}$
#

## Bernacchi set
kc_td_bernacchi = Photosynthesis.KcTDBernacchi(FT)
## CLM set
kc_td_clm       = Photosynthesis.KcTDCLM(FT)

## an array of temperature from 10 to 50 degree C
kc_bernacchi = photo_TD_from_set(kc_td_bernacchi, t_array)
kc_clm       = photo_TD_from_set(kc_td_clm      , t_array)

## plot the TD curves
figure(figsize=(8,6), dpi=100)
tight_layout(true)
plot(t_array, kc_bernacchi, "r-", label="Bernacchi")
plot(t_array, kc_clm      , "g-", label="CLM5"     )
xlabel("Temperature (K)", fontsize=16)
ylabel("\$K_{c}\$ (Pa)" , fontsize=16)
legend(loc="upper left")
gcf()
#----------------------------------------------------------------------------

# ### $K_{o}$
#

## Bernacchi set
ko_td_bernacchi = Photosynthesis.KoTDBernacchi(FT)
## CLM set
ko_td_clm       = Photosynthesis.KoTDCLM(FT)

## an array of temperature from 10 to 50 degree C
ko_bernacchi = photo_TD_from_set(ko_td_bernacchi, t_array)
ko_clm       = photo_TD_from_set(ko_td_clm      , t_array)

## plot the TD curves
figure(figsize=(8,6), dpi=100)
tight_layout(true)
plot(t_array, ko_bernacchi, "r-", label="Bernacchi")
plot(t_array, ko_clm      , "g-", label="CLM5"     )
xlabel("Temperature (K)", fontsize=16)
ylabel("\$K_{o}\$ (Pa)" , fontsize=16)
legend(loc="upper left")
gcf()
#----------------------------------------------------------------------------

# ### $K_{pep}$
#

## Boyd set
kpep_td_boyd = Photosynthesis.KpepTDBoyd(FT)
## CLM set
kpep_td_clm  = Photosynthesis.KpepTDCLM(FT)

## an array of temperature from 10 to 50 degree C
kpep_boyd = photo_TD_from_set(kpep_td_boyd, t_array)
kpep_clm  = photo_TD_from_set(kpep_td_clm , t_array)

## plot the TD curves
figure(figsize=(8,6), dpi=100)
tight_layout(true)
plot(t_array, kpep_boyd, "r-", label="Boyd")
plot(t_array, kpep_clm , "g-", label="CLM5")
xlabel("Temperature (K)" , fontsize=16)
ylabel("\$K_{pep}\$ (Pa)", fontsize=16)
legend(loc="upper left")
gcf()
#----------------------------------------------------------------------------

# ### Respiration
#

## Bernacchi set
re_td_bernacchi = Photosynthesis.RespirationTDBernacchi(FT)
## CLM set
re_td_clm       = Photosynthesis.RespirationTDCLM(FT)

## an array of temperature from 10 to 50 degree C
re_bernacchi = photo_TD_from_val(re_td_bernacchi, FT(1.0), t_array)
re_clm       = photo_TD_from_val(re_td_clm      , FT(1.0), t_array)

## plot the TD curves
figure(figsize=(8,6), dpi=100)
tight_layout(true)
plot(t_array, re_bernacchi, "r-", label="Bernacchi")
plot(t_array, re_clm      , "g-", label="CLM5"     )
xlabel("Temperature (K)"                , fontsize=16)
ylabel("Respiration Rate (μmol m⁻² s⁻¹)", fontsize=16)
legend(loc="upper left")
gcf()
#----------------------------------------------------------------------------

# ### $\Gamma^{*}$
#

## Bernacchi set
Γs_td_bernacchi = Photosynthesis.ΓStarTDBernacchi(FT)
## CLM set
Γs_td_clm       = Photosynthesis.ΓStarTDCLM(FT)

## an array of temperature from 10 to 50 degree C
Γs_bernacchi = photo_TD_from_set(Γs_td_bernacchi, t_array)
Γs_clm       = photo_TD_from_set(Γs_td_clm      , t_array)

## plot the TD curves
figure(figsize=(8,6), dpi=100)
tight_layout(true)
plot(t_array, Γs_bernacchi, "r-", label="Bernacchi")
plot(t_array, Γs_clm      , "g-", label="CLM5"     )
xlabel("Temperature (K)"     , fontsize=16)
ylabel("\$\\Gamma^{*}\$ (Pa)", fontsize=16)
legend(loc="upper left")
gcf()
#----------------------------------------------------------------------------

# ## Summary 
# The way land surface models implement temperature variations in $V_{c,max}$ is still rather variable and there is no clear consensus in what model is best. In addition, gas exchange datasets performed at different temperatures might actually alias different confounding factors into the $V_{c,max}$ determination, so that the actually fitted $V_{c,max}$ strongly depends on the model formulation. E.g. mesophyll conductance is usually ignored in deriving $V_{c,max}$ but might be temperature dependent, which can cause T-dependent errors in $V_{c,max}$. Also, the T-dependencies may change with time, known as acclimation of optimal temperature.
# 
# That said, having the option to use different model forumations (as well as parameter settings) and actually using the ones that are implemented in land surface models directly when fitting leaf level measurements seems to be a prudent way to improve models in the future. More leaf level data might help in the future, especially looking at temperature dependencies, mesophyll conductance, and fluorescence yields. So many chemical and phsyical aspects change with temperature that it is hard to isolate impacting factors. Now, at least, we can run the models directly in the REPL, which is great for prototyping anc validating the model on the individual module level without any effort in terms of additional coding. 

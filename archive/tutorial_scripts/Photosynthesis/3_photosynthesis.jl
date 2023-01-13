#=
# # Leaf Level Photosynthesis
#
# Here, we will go one step further and look at the entire leaf-level response to different environmental conditions, including the impact of leaf diffusive conductance here. See next tutorial for how stomatal conductance responds to the environment.

# ## Leaf Diffusive Conductance and Stomatal Conductance
# Before, we focussed mainly on the demand-driven constraints through Rubisco and RuBP regeneration. Leaf diffusive conductance is highly important as it drives the suppy-side of photosynthesis and interacts with the energy balance and leaf temperature as latent heat fluxes are a major factor in surface cooling.
#
# Before, we have derived net rates of photosynthesis $A_n$, which have to be matched with the supply side through CO₂ diffusion:
#
# $$A_n^{diff} = g_{leaf,C}(C_a-C_c) = g_{leaf,C}\frac{P_a-P_c}{P_{atm}}$$
#
# which can be separated into diffusion from the air to the leaf surface with a boundary layer conductance $g_{b,C}$, diffusion from the surface to the interstitial air-space with stomatal conductance $g_{s,C}$ and diffusion from the  interstitial air-space to the chloroplast with mesophyll conductance $g_{m,C}$ (the letter C stands for CO₂ here, as the diffusion constants vary with species, e.g. H$_2$O stomatal conductance is a factor 1.6 higher than those for CO₂):
#
# $$A_n^{diff} = g_{b,C}(C_a-C_s) = g_{s,C}(C_s-C_i) = g_{m,C}(C_i-C_c) = g_{leaf,C}(C_a-C_c)$$
#
# $$g_{leaf,C} = \left(\frac{1}{g_{b,C}}+\frac{1}{g_{s,C}}+\frac{1}{g_{m,C}}\right)^{-1}$$
#
# The importance of diffusive conductance and its interplay with photosynthesis here is that the supply and demand rates determine $C_c$. A reduction in $C_c$ reduces the demand-limited rates while it increase the diffusion rates (at a given $g$). Most models run internal so-called A-$C_c$ iterations to ensure both rates are balanced and in steady-state ($\partial C_c/\partial t=0$). We implement this option but also opt to run stomatal conductance prognostically, as stomata don't open and close instantanously but have a response time of around 15 minutes.
#
# Below, we show simple examples of how the environmental conditions and the leaf diffusive conductance impact leaf photosynthesis.
#

## load packages
using PyPlot
using Land.Photosynthesis

## set the default floating type
FT     = Float32;

## define default photosynthesis parameter sets
c3_set = C3CLM(FT);
c4_set = C4CLM(FT);

## define leaf photosynthetic parameters and environmental conditions
leaf_3 = LeafPhotoContainer{FT}(APAR=1200, Vcmax25=60, Jmax25=120, Vpmax25=80, Rd25=1);
leaf_4 = LeafPhotoContainer{FT}(APAR=1200, Vcmax25=60, Jmax25=120, Vpmax25=80, Rd25=1);
envir  = EnvironmentConditions{FT}();
#----------------------------------------------------------------------------

# ## $A$ ~ $P_c$

## set the internal CO₂ to 30 Pa and temperature to 298.15 K and PAR to 1500
leaf_3.p_i=30; leaf_3.T=298.15; leaf_3.APAR=1500;
leaf_4.p_i=30; leaf_4.T=298.15; leaf_4.APAR=1500;

## set an array of internal CO₂
p_array = collect(FT, 5:200);
ac_3 = similar(p_array); aj_3 = similar(p_array); ap_3 = similar(p_array); ag_3 = similar(p_array); an_3 = similar(p_array);
ac_4 = similar(p_array); aj_4 = similar(p_array); ap_4 = similar(p_array); ag_4 = similar(p_array); an_4 = similar(p_array);

## update the temperature and PAR dependent photosynthesis
Photosynthesis.photosynthesis_TD!(c3_set, leaf_3, envir);
Photosynthesis.photosynthesis_TD!(c4_set, leaf_4, envir);


## update the CO₂ dependent photosynthesis
for i in eachindex(p_array);
    p_i = p_array[i];
    photosynthesis_PD!(c3_set, leaf_3, p_i);
    photosynthesis_PD!(c4_set, leaf_4, p_i);

    ac_3[i], aj_3[i], ap_3[i], ag_3[i], an_3[i] = leaf_3.Ac, leaf_3.Aj, leaf_3.Ap, leaf_3.Ag, leaf_3.An;
    ac_4[i], aj_4[i], ap_4[i], ag_4[i], an_4[i] = leaf_4.Ac, leaf_4.Aj, leaf_4.Ap, leaf_4.Ag, leaf_4.An;
end

## plot the CO₂ responses for C3 and C4 photosynthesis
figure(figsize=(10,6), dpi=100);
tight_layout(true);
subplot(1,2,1);
plot(p_array, ac_3, "r:", lw=2, label="Ac");
plot(p_array, aj_3, "b:", lw=2, label="Aj");
plot(p_array, ap_3, "g:", lw=2, label="Ap");
plot(p_array, ag_3, "k-", lw=1, label="Ag");
plot(p_array, an_3, "k:", lw=1, label="An");
xlabel("Leaf Internal CO₂ (Pa)"            , fontsize=16);
ylabel("Photosynthetic Rate (μmol m⁻² s⁻¹)", fontsize=16);
legend(loc="lower right")
subplot(1,2,2);
plot(p_array, ac_4, "r:", lw=2, label="Ac");
plot(p_array, aj_4, "b:", lw=2, label="Aj");
plot(p_array, ap_4, "g:", lw=2, label="Ap");
plot(p_array, ag_4, "k-", lw=1, label="Ag");
plot(p_array, an_4, "k:", lw=1, label="An");
xlabel("Leaf Internal CO₂ (Pa)"            , fontsize=16);
legend(loc="lower right")
gcf()
#----------------------------------------------------------------------------

# ## $A$ ~ $PAR$

## set the internal CO₂ to 30 Pa and temperature to 298.15 K and PAR to 1500
leaf_3.p_i=30; leaf_3.T=298.15; leaf_3.APAR=1500;
leaf_4.p_i=30; leaf_4.T=298.15; leaf_4.APAR=1500;

## set an array of APAR
par_array = collect(FT, 10:10:2000);
ac_3 = similar(par_array); aj_3 = similar(par_array); ap_3 = similar(par_array); ag_3 = similar(par_array); an_3 = similar(par_array);
ac_4 = similar(par_array); aj_4 = similar(par_array); ap_4 = similar(par_array); ag_4 = similar(par_array); an_4 = similar(par_array);

for i in eachindex(par_array);
    leaf_3.APAR = par_array[i];
    leaf_4.APAR = par_array[i];

    ## update the temperature dependent photosynthesis
    Photosynthesis.photosynthesis_TD!(c3_set, leaf_3, envir);
    Photosynthesis.photosynthesis_TD!(c4_set, leaf_4, envir);

    ## update the CO₂ dependent photosynthesis
    photosynthesis_PD!(c3_set, leaf_3, leaf_3.p_i);
    photosynthesis_PD!(c4_set, leaf_4, leaf_4.p_i);

    ac_3[i], aj_3[i], ap_3[i], ag_3[i], an_3[i] = leaf_3.Ac, leaf_3.Aj, leaf_3.Ap, leaf_3.Ag, leaf_3.An;
    ac_4[i], aj_4[i], ap_4[i], ag_4[i], an_4[i] = leaf_4.Ac, leaf_4.Aj, leaf_4.Ap, leaf_4.Ag, leaf_4.An;
end

## plot the CO₂ responses for C3 and C4 photosynthesis
figure(figsize=(10,6), dpi=100);
tight_layout(true);
subplot(1,2,1);
plot(par_array, ac_3, "r:", lw=2, label="Ac");
plot(par_array, aj_3, "b:", lw=2, label="Aj");
plot(par_array, ap_3, "g:", lw=2, label="Ap");
plot(par_array, ag_3, "k-", lw=1, label="Ag");
plot(par_array, an_3, "k:", lw=1, label="An");
xlabel("Absorbed PAR (μmol m⁻² s⁻¹)"       , fontsize=16);
ylabel("Photosynthetic Rate (μmol m⁻² s⁻¹)", fontsize=16);
legend(loc="lower right")
subplot(1,2,2);
plot(par_array, ac_4, "r:", lw=2, label="Ac");
plot(par_array, aj_4, "b:", lw=2, label="Aj");
plot(par_array, ap_4, "g:", lw=2, label="Ap");
plot(par_array, ag_4, "k-", lw=1, label="Ag");
plot(par_array, an_4, "k:", lw=1, label="An");
xlabel("Absorbed PAR (μmol m⁻² s⁻¹)"       , fontsize=16);
legend(loc="lower right")
gcf()
#----------------------------------------------------------------------------

# ## $A$ ~ $T$

## set the internal CO₂ to 30 Pa and temperature to 298.15 K and PAR to 1500
leaf_3.p_i=30; leaf_3.T=298.15; leaf_3.APAR=1500;
leaf_4.p_i=30; leaf_4.T=298.15; leaf_4.APAR=1500;

## set an array of APAR
t_array = collect(FT, 10:50) .+ FT(273.15);
ac_3 = similar(t_array); aj_3 = similar(t_array); ap_3 = similar(t_array); ag_3 = similar(t_array); an_3 = similar(t_array);
ac_4 = similar(t_array); aj_4 = similar(t_array); ap_4 = similar(t_array); ag_4 = similar(t_array); an_4 = similar(t_array);

for i in eachindex(t_array);
    leaf_3.T = t_array[i];
    leaf_4.T = t_array[i];

    ## update the temperature dependent photosynthesis
    Photosynthesis.photosynthesis_TD!(c3_set, leaf_3, envir);
    Photosynthesis.photosynthesis_TD!(c4_set, leaf_4, envir);

    ## update the CO₂ dependent photosynthesis
    photosynthesis_PD!(c3_set, leaf_3, leaf_3.p_i);
    photosynthesis_PD!(c4_set, leaf_4, leaf_4.p_i);

    ac_3[i], aj_3[i], ap_3[i], ag_3[i], an_3[i] = leaf_3.Ac, leaf_3.Aj, leaf_3.Ap, leaf_3.Ag, leaf_3.An;
    ac_4[i], aj_4[i], ap_4[i], ag_4[i], an_4[i] = leaf_4.Ac, leaf_4.Aj, leaf_4.Ap, leaf_4.Ag, leaf_4.An;
end

## plot the CO₂ responses for C3 and C4 photosynthesis
figure(figsize=(10,6), dpi=100);
tight_layout(true);
subplot(1,2,1);
plot(t_array, ac_3, "r:", lw=2, label="Ac");
plot(t_array, aj_3, "b:", lw=2, label="Aj");
plot(t_array, ap_3, "g:", lw=2, label="Ap");
plot(t_array, ag_3, "k-", lw=1, label="Ag");
plot(t_array, an_3, "k:", lw=1, label="An");
xlabel("Temperature (K)"                   , fontsize=16);
ylabel("Photosynthetic Rate (μmol m⁻² s⁻¹)", fontsize=16);
legend(loc="upper left")
subplot(1,2,2);
plot(t_array, ac_4, "r:", lw=2, label="Ac");
plot(t_array, aj_4, "b:", lw=2, label="Aj");
plot(t_array, ap_4, "g:", lw=2, label="Ap");
plot(t_array, ag_4, "k-", lw=1, label="Ag");
plot(t_array, an_4, "k:", lw=1, label="An");
xlabel("Temperature (K)"                   , fontsize=16);
legend(loc="upper left")
gcf()
#----------------------------------------------------------------------------
=#

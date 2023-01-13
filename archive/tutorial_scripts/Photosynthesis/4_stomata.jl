#=
# # Stomatal conductance
#
# Here, we will go one step further and look at the entire leaf-level response to different environmental conditions, including the impact of stomatal responses to the environment here.
#
# ## Empirical Stomatal Conductance models
#
# Currently, we can choose between the widely used Ball-Berry model[^1] or the Medlyn approach[^2].
#
# Ball and Berry derived the following empirical formulation based on leaf level data (RH is relative humidity):
#
# $$g_{s,w} = g_0 + g_1  \frac{A_n \cdot RH}{C_s}$$
#
# Medlyn derived the following equations based on stomatal optimization theory but in a similar form as the BB model:
#
# $$g_{s,w} = g_0 + \left(1 + \frac{g_1}{\sqrt{VPD}}\right) \frac{A_n}{C_s}$$
#
# Both models are proportional to $A_n$ and inversely proportional to $C_s$, with the main difference in the dependence on either relative humidity in the Ball-Berry model vs. vapor pressure deficit (VPD) in Medlyn et al. In both cases, $g_0$ is the residual conductivity even if stomata are fully closed and $g_1$ is related to the marginal water cost of plant carbon gain. Importantly, $g_1$ can't be inter-changed between the formulations, underlining again that parameters have to be optimized with respect to the model that is eventually being used.
#
# ## Stomatal Optimization Theories
#
# The empirical formulations only hold in well-watered soil and our main goal is to implement stomatal optimization models to take the entire soil-plant-atmosphere continuum into account[^3]. Here, we will just use the empirical models and steady-state photosynthesis to show the underlying principles.
#
# [^1]: Ball, J.T., Woodrow, I.E. and Berry, J.A., 1987. A model predicting stomatal conductance and its contribution to the control of photosynthesis under different environmental conditions. In Progress in photosynthesis research (pp. 221-224). Springer, Dordrecht.
#
# [^2]: Medlyn, B.E., Duursma, R.A., Eamus, D., Ellsworth, D.S., Prentice, I.C., Barton, C.V., Crous, K.Y., De Angelis, P., Freeman, M. and Wingate, L., 2011. Reconciling the optimal and empirical approaches to modelling stomatal conductance. Global Change Biology, 17(6), pp.2134-2144.
#
# [^3]: Wang, Y., Sperry, J.S., Anderegg, W.R., Venturas, M.D. and Trugman, A.T., 2020. A theoretical and empirical assessment of stomatal optimization modeling. New Phytologist.
#

## add usual tools
using BenchmarkTools
using PyPlot

## load photosynthesis modules
using Land.WaterPhysics
using Land.Photosynthesis
using Land.Plant

## specify floating type
FT = Float32;
#----------------------------------------------------------------------------

## define default photosynthesis parameter sets
c3_set = C3CLM(FT);
c4_set = C4CLM(FT);

## define leaf photosynthetic parameters and environmental conditions
leaf_3 = LeafPhotoContainer{FT}(APAR=1200, Vcmax25=60, Jmax25=120, Vpmax25=80, Rd25=1);
leaf_4 = LeafPhotoContainer{FT}(APAR=1200, Vcmax25=60, Jmax25=120, Vpmax25=80, Rd25=1);
envir  = EnvironmentConditions{FT}();

## create empirical stomatal scheme
sto_c3 = ESMBallBerry{FT}(g1 = 20)
sto_c4 = ESMBallBerry{FT}(g1 = 4)
#----------------------------------------------------------------------------

# ## Stomatal response to CO₂

## set the internal CO₂ to 30 Pa and temperature to 298.15 K and PAR to 1500
leaf_3.p_i=30; leaf_3.T=298.15; leaf_3.APAR=1500;
leaf_4.p_i=30; leaf_4.T=298.15; leaf_4.APAR=1500;

## set an array of internal CO₂
p_array = collect(FT, 5:200);
an_3 = similar(p_array); gs_3 = similar(p_array); pi_3 = similar(p_array);
an_4 = similar(p_array); gs_4 = similar(p_array); pi_4 = similar(p_array);

for i in eachindex(p_array)
    envir.p_a = p_array[i];
    get_empirical_gsw_pi(c3_set, leaf_3, envir, sto_c3);
    get_empirical_gsw_pi(c4_set, leaf_4, envir, sto_c4);
    an_3[i], gs_3[i], pi_3[i] = leaf_3.An, leaf_3.g_sw, leaf_3.p_i
    an_4[i], gs_4[i], pi_4[i] = leaf_4.An, leaf_4.g_sw, leaf_4.p_i
end

## plot the results
figure(figsize=(12,3.5), dpi=100)
tight_layout(true)
subplot(1,3,1)
plot(p_array, an_3, label="C3")
plot(p_array, an_4, label="C4")
xlabel("Leaf internal CO₂ (Pa)", fontsize=12)
ylabel("Anet (μmol m⁻² s⁻¹)"   , fontsize=12)
legend()
subplot(1,3,2)
plot(p_array, gs_3, label="C3")
plot(p_array, gs_4, label="C4")
xlabel("Leaf internal CO₂ (Pa)"            , fontsize=12)
ylabel("Stomatal conductance (mol m⁻² s⁻¹)", fontsize=12)
legend()
subplot(1,3,3)
plot(p_array, pi_3 ./ p_array, label="C3")
plot(p_array, pi_4 ./ p_array, label="C4")
xlabel("Leaf internal CO₂ (Pa)", fontsize=12)
ylabel("Cc/Ca"                 , fontsize=12)
legend()
gcf()
#----------------------------------------------------------------------------

# ## Stomatal response to PAR

## set the internal CO₂ to 30 Pa and temperature to 298.15 K and PAR to 1500
leaf_3.p_i=30; leaf_3.T=298.15; leaf_3.APAR=1500;
leaf_4.p_i=30; leaf_4.T=298.15; leaf_4.APAR=1500;
envir  = EnvironmentConditions{FT}();

## set an array of internal CO₂
par_array = collect(FT, 0:10:2000);
an_3 = similar(par_array); gs_3 = similar(par_array); pi_3 = similar(par_array);
an_4 = similar(par_array); gs_4 = similar(par_array); pi_4 = similar(par_array);

for i in eachindex(par_array)
    leaf_3.APAR = par_array[i];
    leaf_4.APAR = par_array[i];
    get_empirical_gsw_pi(c3_set, leaf_3, envir, sto_c3);
    get_empirical_gsw_pi(c4_set, leaf_4, envir, sto_c4);
    an_3[i], gs_3[i], pi_3[i] = leaf_3.An, leaf_3.g_sw, leaf_3.p_i
    an_4[i], gs_4[i], pi_4[i] = leaf_4.An, leaf_4.g_sw, leaf_4.p_i
end

## plot the results
figure(figsize=(12,3.5), dpi=100)
tight_layout(true)
subplot(1,3,1)
plot(par_array, an_3, label="C3")
plot(par_array, an_4, label="C4")
xlabel("PAR (μmol m⁻² s⁻¹)" , fontsize=12)
ylabel("Anet (μmol m⁻² s⁻¹)", fontsize=12)
legend()
subplot(1,3,2)
plot(par_array, gs_3, label="C3")
plot(par_array, gs_4, label="C4")
xlabel("PAR (μmol m⁻² s⁻¹)"                , fontsize=12)
ylabel("Stomatal conductance (mol m⁻² s⁻¹)", fontsize=12)
legend()
subplot(1,3,3)
plot(par_array, pi_3 ./ envir.p_a, label="C3")
plot(par_array, pi_4 ./ envir.p_a, label="C4")
xlabel("PAR (μmol m⁻² s⁻¹)", fontsize=12)
ylabel("Cc/Ca"             , fontsize=12)
legend()
gcf()
#----------------------------------------------------------------------------

# ## Stomatal response to Temperature

## set the internal CO₂ to 30 Pa and temperature to 298.15 K and PAR to 1500
leaf_3.p_i=30; leaf_3.T=298.15; leaf_3.APAR=1500;
leaf_4.p_i=30; leaf_4.T=298.15; leaf_4.APAR=1500;
envir  = EnvironmentConditions{FT}();

## set an array of internal CO₂
t_array = collect(FT, 280:320);
an_3 = similar(t_array); gs_3 = similar(t_array); pi_3 = similar(t_array);
an_4 = similar(t_array); gs_4 = similar(t_array); pi_4 = similar(t_array);

for i in eachindex(t_array)
    leaf_3.T = t_array[i];
    leaf_4.T = t_array[i];
    envir.p_H₂O = saturation_vapor_pressure(t_array[i]) * 0.5;
    get_empirical_gsw_pi(c3_set, leaf_3, envir, sto_c3);
    get_empirical_gsw_pi(c4_set, leaf_4, envir, sto_c4);
    an_3[i], gs_3[i], pi_3[i] = leaf_3.An, leaf_3.g_sw, leaf_3.p_i
    an_4[i], gs_4[i], pi_4[i] = leaf_4.An, leaf_4.g_sw, leaf_4.p_i
end

## plot the results
t_array .-= FT(273.15)

figure(figsize=(12,3.5), dpi=100)
tight_layout(true)
subplot(1,3,1)
plot(t_array, an_3, label="C3")
plot(t_array, an_4, label="C4")
xlabel("Temperature (K)"    , fontsize=12)
ylabel("Anet (μmol m⁻² s⁻¹)", fontsize=12)
legend()
subplot(1,3,2)
plot(t_array, gs_3, label="C3")
plot(t_array, gs_4, label="C4")
xlabel("Temperature (K)"                   , fontsize=12)
ylabel("Stomatal conductance (mol m⁻² s⁻¹)", fontsize=12)
legend()
subplot(1,3,3)
plot(t_array, pi_3 ./ envir.p_a, label="C3")
plot(t_array, pi_4 ./ envir.p_a, label="C4")
xlabel("Temperature (K)", fontsize=12)
ylabel("Cc/Ca"          , fontsize=12)
legend()
gcf()
#----------------------------------------------------------------------------




##@btime get_empirical_gsw_pi(c3_set, leaf_3, envir, sto_c3);
##@btime get_empirical_gsw_pi(c4_set, leaf_4, envir, sto_c4);
=#

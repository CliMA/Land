
# # Leaf Photosynthesis Basics
# This tutorial will walk you through the most basic aspects of how we implement Photosynthesis at the leaf level. Most of the concepts are described in the literature, with the first quantitative approach to modeling photosynthesis given in Farquhar, von Caemmerer and Berry[^1] in their seminal 1980 paper for C3 photosynthesis. C4 photosynthesis is largely based on Collatz et al[^2] but we approach the photosynthesis modeling with a rather generic approach that facilitates the application of different photosynthesis modeling approaches without changing the core code base. A good overview on the entire process of photosynthesis and different parameterizations can be found in Bonan[^3]. 
# 
# At the core of both C3 and C4 photosynthesis is an enzyme catalyzed reaction of ribulose-1,5-bisphosphate (RuBP) with 
# CO$_2$, yielding two 3-carbon compounds (phosphoglycerate (PGA)) as the initial products of photosynthesis. The enzyme RuBP carboxylase/oxygenase (Rubisco) catalyzes this reaction. With the regeneration of the substrate RuBP through the light reactions (using produced ATP and NADPH), the core cycle of photosynthesis is formed. This cycle was discovered in 1950 by Melvin Calvin, James Bassham, and Andrew Benson at the University of California, Berkeley[^4][^5] by using the radioactive isotope carbon-14.
# 
# An oxygenation step of RuBP releases half a CO$_2$. This so called photorespiration process results in inefficiencies in C3 photosynthesis, as both CO$_2$ and O$_2$ compete at the Rubisco site. The overall photosynthetic rate of the enzyme-catalyzed turnover rates at the Rubisco site thus determine the overall photosynthesis
# 
# $$A_n = V_c - 0.5V_o - R_d\,,$$
# 
# with $V_c$ being the carboxylation rate, $V_o$ the oxygenation rate and $R_d$ the mitochondrial respiration. Both rates follow Michaelis-Menten kinetics, accounting for the competing substrate effects:
# 
# $$V_c = \frac{V_{c,max}C_c}{C_c+K_c(1+O_c/K_o)}$$
# 
# $$V_o = \frac{V_{o,max}O_c}{O_c+K_o(1+C_c/K_c)}\,$$
# 
# with $K_c$ and $K_o$ being the Michaelis Menten constants for CO$_2$ and O$_2$, $C_c$ and $O_c$ the partial pressures of CO$_2$ and O$_2$ at the Rubisco site.
# 
# The ratio of oxygenation to carboxylation rates is
# 
# $$\phi = \frac{V_o}{V_c} = \frac{V_{o,max}K_c}{V_{c,max}K_o}\frac{O_c}{C_c}\,$$ 
# 
# which yields the CO$_2$ compensation point $\Gamma_\star$:
# 
# $$\Gamma_\star = 0.5\frac{V_{o,max}K_c}{V_{c,max}K_o}O_c\,$$ 
# 
# with is the internal CO$_2$ partial pressure at which oxygenation and carboxylation cancel each other out in terms of CO$_2$ consumption and production (neglecting $R_d$).
# 
# 
# 
# 
# [^1]: Farquhar, G.D., von Caemmerer, S.V. and Berry, J.A., 1980. A biochemical model of photosynthetic CO$_2$ assimilation in leaves of C3 species. Planta, 149(1), pp.78-90.
# 
# [^2]: Collatz, G.J., Ribas-Carbo, M. and Berry, J.A., 1992. Coupled photosynthesis-stomatal conductance model for leaves of C4 plants. Functional Plant Biology, 19(5), pp.519-538.
# 
# [^3]: Bonan, G., 2019. Climate change and terrestrial ecosystem modeling. Cambridge University Press.
# 
# [^4]: Calvin, Melvin, and Andrew Alm Benson. "The path of carbon in photosynthesis IV: the identity and sequence of the intermediates in sucrose synthesis." Science 109, no. 2824 (1949): 140-142.
# 
# [^5]: Benson, A.A., Bassham, J.A., Calvin, M., Goodale, T.C., Haas, V.A. and Stepka, W., 1950. The path of carbon in photosynthesis. v. paper chromatography and radioautography of the products1. Journal of the American Chemical Society, 72(4), pp.1710-1718.
# 

# ## Rate Limiting Steps for photosynthesis
# 
# ### Rubisco-limited rates:
# The net photosynthetic rate limited by Rubisco when RuBP re-generation is not constraining can thus be described as 
# 
# $$A_n = \left(1-\frac{\Gamma_\star}{C_c}\right)V_c-R_d\,$$
# 
# which equals (now denoting $A_c$ as the Rubisco limited rate):
# 
# $$A_c = \frac{(C_c-\Gamma_\star) V_{c,max}}{C_c+K_c(1+O_c/K_o)}\,$$
# 
# which is implemented in our routines "rubisco_limited_rate!".
# 
# ### RuBP-regeration limited rates (light-limited):
# 
# RuBP is regenerated via the ligh reaction, which generates ATP and NADP to power this part of the Calvin-Benson-Bassham cycle. Without going into details of the NADPH or ATP requirements for regeneration, the rate of RuBP limited photosynthesis through light-powered electron transport $J$ (µmol/m$^2$/s) is given as 
# 
# $$A_j = \underbrace{\frac{C_c-\Gamma_\star}{Cc}}_{\text{loss in photorespiration}}\underbrace{\frac{J\,C_c}{4C_c+8\Gamma_\star}}_{\text{RuBP regeneration}} = \frac{J(C_c-\Gamma_\star)}{4C_c+8\Gamma_\star}\,$$
# 
# One can already see what the key difference between C3 and C4 photosynthesis is based on these set of equation for rubisco turnover and RuBP-regeneration limited rates. In C4 photosynthesis, $C_c$ is typically only about 70% of the ambient CO$_2$ concentration as CO$_2$ has to diffuse through stomata and the mesophyll. For C4, a carbon accumulation mechanism uses a four-carbon organic acid (hence C4) to transport CO$_2$ from the mesophyll cells to the bundle sheath cells, where the reaction with Rubisco takes place. This leads to much higher $C_c$ for C4 plants ($>>$ than ambient air CO$_2$), which is often simplified for $A_c$ as the limit of C3 equation with $\lim_{C_c \to \infty}$ Ac(C3):
# 
# $$A_c^{C4} = V_{c,max}$$
# 
# Similarly, the RuBP-regeneration limited rate simplifies to:
# 
# $$A_j^{C4} = \alpha J$$
# 
# with a slightly lower efficiency $\alpha$ for C4 photosynthesis compared to C3 (as the carbon accumulation mechanism also consumes ATP). In reality, the situation is somewhat more complex and there are different approximations for Rubisco and RuBP limited rate constants. A comprehensive overview is described in von Caemmerer [^6]. Again, our goal is to provide a flexible framework for photosynthesis modeling, so different implementations of rate limiting steps $A_c,A_j$ can be used in a modular framework (achieved through code abstraction and multiple dispatch in Julia).
# 
# ### Electron Transport Rate $J$
# 
# The rate of electron transport is driven by absorbed photosynthetically absorbed radiation (APAR, $\mu mol/m^2/s$) by Photosystem II (PSII). With our leaf level optical model, we already compute the efficiency of absorbtion by leaves depending on pigment contents. We define $\varphi_{PSII}$ as the quantum yield of photosystem II (maxima about 0.83) and $f_{PSII}$ the fraction of light used for PSII (PSII/(PSI+PSII)), typically assume to be 0.5.
# 
# $$J_{PSII} = f_{PSII} \varphi_{PSII} APAR$$
# 
# In most models, a maximum electron transport $J_{max}$ rate is assumed and the actual electron trapsort rate $J_a$ is defined as the lower root of the quadratic expression
# 
# $$\Theta_j J_a^2 -(J_{PSII}+J_{max})+ J_{PSII}J_{max} = 0$$
# 
# where $\Theta_j$ is a curvature parameter to assure a smooth transition. 
# 
# ### Product limited rates:
# 
# Typically, a thrid limitation of photosynthesis is being used as well, which we denote as product limited rate here, even though the processes for C3 and C4 plants differ. Typically, the export of the products of photosynthesis (triose phosphates) in the synthesis of sugars can be rate limiting, which is often parameterized by $V_{c,max}$:
# 
# $$A_p = a V_{c,max}\,,$$
# 
# where our standard definition for C3 used a=0.5, which rarely limits photosynthesis compared to $A_c$ and $A_j$.
# 
# For C4 plants, we use the PEP-carboxylase CO$_2$ concentration mechanism into the bundle sheath cell as product limited rate step, using Michaelis-Menten kinetics:
# 
# $$A_p^{C4} = \frac{V_{p,max}C_c}{C_c + K_p}$$
# 
# with corresponding Michaelis Menten constant $K_p$ and maximal rates $V_{p,max}$.
# 
# 
# ### Total rate:
# 
# The gross photosynthetic rate $A_g$ can then be define as the mininum of all possible limitations:
# 
# $$A_g = min(Ac,Aj,Ap)$$
# 
# or, alternatively, using quadratic equations as for $J$, which provides smoother transitions and co-limitation to some degree (both are options in our setup with a user defined curvature parameter $\Theta$.
# 
# Net leaf photosynthesis includes mitochondrial respiration as well:
# 
# $$A_n = A_g - R_d$$
# 
# [^6]: Von Caemmerer, S., 2000. Biochemical models of leaf photosynthesis. Csiro publishing.

# ## Simple example:
# We just want to lay out a first simple example as to how our model setup work:  

## Loading the Photosynthesis model:
using Revise
using Land.Photosynthesis
## Defining our Field Type (we can easily switch between Double and Float precision this way)
const FT = Float32;
#----------------------------------------------------------------------------

## Create a standard leaf with defualt parameters
leaf = leaf_params{FT}();

## Create a standard meteo structure:
met = meteo{FT}();
#----------------------------------------------------------------------------

# ----
# How to use the documentation, what do we know about leaf_params, which stores most physiologically relevant parameters. 

?leaf_params
#----------------------------------------------------------------------------

?meteo
#----------------------------------------------------------------------------

# ### Defining the model setup
# The most important step is to define which submodules to use. There might be different implementations for Fluorescence, Photosynthesis (C3,C4,etc), respiration, stomatal conductance (e.g. Ball-Berry, Medlyn), the T-dependence of J$_{max}$, V$_{c,max}$ and Michaelis Menten constants as well as leaf boundary layer resistance (setting it to 0 here to mimic well vented laboratory leaf level measurements) and colimitation.   

##Here, we just have to worry abou the photosynthesis module, which we set here:
mod_photo = C3FvCBPhoto()

## All modules here:
mods = Photosynthesis.PhotoMods(
    fluorescence    = FlexasTolBerryFluorescence{FT}(),
    photosynthesis  = C3FvCBPhoto(),
    respiration     = RespirationCLM{FT}(),
    stomatal        = BallBerryStomata{FT}(g1=8),
    Jmax            = JmaxCLM{FT}(),
    Vmax            = VcmaxCLM{FT}(),
    MichaelisMenten = MM_CLM{FT}(),
    BoundaryLayer   = FixedBoundaryResistance{FT}(ra=0),
    colimitation = CurvedColimit{FT}(0.99));
#----------------------------------------------------------------------------

## Set APAR to 250 $\mu mol/m^2/s$
leaf.APAR = 250;
## Set temperature to 290K
leaf.T = 290;
## Applying the T-correction for all rate constants:
Photosynthesis.set_leaf_temperature!(mods, leaf)
@show leaf.Vcmax
@show leaf.Jmax
## Specify Cc directly here in ppm (will be converted to Pa internally)
leaf.Cc = 350;
#----------------------------------------------------------------------------

Aj = Photosynthesis.light_limited_rate!(mod_photo, leaf, met, leaf.APAR)
#----------------------------------------------------------------------------

Ac = Photosynthesis.rubisco_limited_rate!(mod_photo, leaf, met)
#----------------------------------------------------------------------------

Ap = Photosynthesis.product_limited_rate!(mod_photo, leaf, met)
#----------------------------------------------------------------------------

# ----
# ### Summary
# This was the simplest example to compute A$_c$, A$_j$ and A$_p$ with a specified C$_c$ and APAR. In reality, C$_c$ is defined by the interplay of photosynthetic demand and the supply through diffusion through stomata and the mesophyll. However, you can easily exchange the photosynthesis model above with C4FvCBPhoto(). You can also play around with leaf temperature T and APAR as well as $C_c$. 

?C3FvCBPhoto
#----------------------------------------------------------------------------

?C4CollatzPhoto
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------

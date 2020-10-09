# Photosynthesis

```@meta
CurrentModule = Land.Photosynthesis
```




## Leaf and Environment Structures
To model photosynthesis more efficiently, we use a container ([`Leaf`](@ref)
    struct) to store the photosynthesis-related information. For example, many
    of the physiological parameters are temperature-dependent, but these
    temperature-dependent values only need to be updated when leaf temperature
    changes. Therefore, use of the container significantly reduces the time
    required when programing leaf gas exchange prognostically. The
    [`Leaf`](@ref) struct has the following fields:

```@docs
Leaf
```

Also, environmental conditions are required to compute photosynthetic rate, and
    these conditions are stored in [`AirLayer`](@ref) struct. An
    [`AirLayer`](@ref) struct further allows for more conveniently modeling
    photosynthesis the vertical CO₂ and H₂O gradients in the canopy. The
    [`AirLayer`](@ref) structs has the following fields:

```@docs
AirLayer
```

See exmaples below for how to create the structs
```julia
using Photosynthesis

FT = Float32;
leaf = Leaf{FT}();
envir = AirLayer{FT}();
```




## Temperature Dependency Structs
The temperature-dependent (TD) photosynthetic parameters include
- ``J_\text{max}`` Maximal electron transport rate
- ``K_\text{c}`` Michaelis constant for CO₂
- ``K_\text{m}`` Michaelis-Menten's coefficient
- ``K_\text{o}`` Michaelis constant for O₂
- ``K_\text{pep}`` Michaelis constant for PEP carboxylation
- ``R_\text{d}`` Dark respiration
- ``V_\text{cmax}`` Maximal RuBP carboxylation rate
- ``V_\text{omax}`` Maximal RuBP oxygenation rate
- ``V_\text{pmax}`` Maximal PEP carboxylation rate
- ``Γ^{*}`` CO₂ compensation point with the absence of dark respiration

There are two typical types of temperature dependencies using the classic
    Arrhenius equation. We define the three types as [`ArrheniusTD`](@ref),
    [`ArrheniusPeakTD`](@ref), and [`Q10TD`](@ref) subject to
    [`AbstractTDParameterSet`](@ref) type:

```@docs
AbstractTDParameterSet
ArrheniusTD
ArrheniusPeakTD
Q10TD
```

There are many published parameter sets for the various temperature
    dependencies, and to ease the modeling we predefined most of the structs:
```@docs
JmaxTDBernacchi
JmaxTDCLM
JmaxTDLeuning
KcTDBernacchi
KcTDCLM
KoTDBernacchi
KoTDCLM
KpepTDBoyd
KpepTDCLM
RespirationTDBernacchi
RespirationTDCLM
VcmaxTDBernacchi
VcmaxTDCLM
VcmaxTDLeuning
VomaxTDBernacchi
VpmaxTDBoyd
ΓStarTDBernacchi
ΓStarTDCLM
```

The TDs can be easily created using commands like

```julia
using Photosynthesis

FT = Float32;
_td_1 = JmaxTDBernacchi(FT);
_td_2 = VcmaxTDCLM(FT);
```

However, be aware that these pre-defined TD structs are not mutable, to create
    customized TD struct, code like this will be useful

```julia
using Photosynthesis

FT = Float32;
_td_1 = ArrheniusTD{FT}(1, 10000, 30);
_td_1 = ArrheniusPeakTD{FT}(1, 10000, 30, 1);
```

To further simplify the use of Photosynthesis module, we provide a few
    collections/structs of temperature dependencies as well as other parameter
    sets like [`FluoParaSet`](@ref). The structs are catergorized to
    [`C3ParaSet`](@ref) and [`C4ParaSet`](@ref) subject to an
    [`AbstractPhotoModelParaSet`](@ref) type, and the structs are meant for
    modeling C3 photosynthesis and C4 photosynthesis, respectively.

```@docs
AbstractPhotoModelParaSet
C3ParaSet
C4ParaSet
```

Again, to guarantee a quick start, we provided a few pre-defined parameter
    sets:

```@docs
C3Bernacchi
C3CLM
C4CLM
```

Examples:
```julia
using Photosynthesis

FT = Float32;
set_b = C3Bernacchi(FT);
set_3 = C3CLM(FT);
set_4 = C4CLM(FT);
```

Note it here that the [`C3ParaSet`](@ref) and [`C4ParaSet`](@ref) structs are
    mutable, and the fields can be changed to another non-mutable TD struct.
    We'd like to mention that in some cases, leaf respiration rate is not
    measured, and in this case, the dark respiration rate will be computed from
    $V_\text{cmax}$ using a multiplier

```@docs
VtoRCollatz
VtoRDefault
```




## Temperature Dependency
As mentioned above, temperature corrections only need to be done once per
    temperature change, and storing the temperature corrected values will
    significantly boost the code speed. Here we provide a few functions to
    change the stored values. First of all, all the temperature corrections are
    made with [`temperature_correction`](@ref):

```@docs
temperature_correction
```

Second, depending on which physiological parameter to correct, some corrections
    use the `VAL_25` field in the [`ArrheniusTD`](@ref), like $K_\text{c}$,
    $K_\text{o}$, and $K_\text{pep}$:

```@docs
photo_TD_from_set
```

Some corrections use the reference values from the [`Leaf`](@ref) struct, like
    $V_\text{cmax}$ and $J_\text{max}$:

```@docs
photo_TD_from_val
```

The functions to make temperature corrections to each individual variables are

```@docs
leaf_jmax!
leaf_kc!
leaf_km!
leaf_ko!
leaf_kpep!
leaf_rd!
leaf_vcmax!
leaf_vpmax!
leaf_Γstar!
```

Again to ease the coding, we provide a function to run all the temperature
    dependencies:

```@docs
leaf_temperature_dependence!
```

Note it here that function [`leaf_temperature_dependence!`](@ref) updates
    saturated vapor pressure from leaf temperature as well.

Example:

```julia
using Photosynthesis

FT = Float32;
leaf = Leaf{FT}();
envir = AirLayer{FT}();
set_3 = C3CLM(FT);

leaf_temperature_dependence!(c3_set, leaf, envir);
leaf_temperature_dependence!(c3_set, leaf, envir, FT(300));
```




## RubisCO-limited Photosynthesis
By default, Photosynthesis module computes gross photosynthetic rate as the
    minimal of the three:

- ``A_\text{c}`` RubisCO-limited photosynthetic rate
- ``A_\text{j}`` Light-limited photosynthetic rate
- ``A_\text{p}`` Product-limited photosynthetic rate

If leaf internal CO₂ is known, $A_\text{c}$ (gross rate) can be computed using

```@docs
rubisco_limited_rate!
```

If total leaf diffusive conductance to CO₂ is known, $A_\text{c}$ can be
    computed analytically by solving the quadratic function using

```@docs
lower_quadratic
```

The function to analytically compute $A_\text{c}$ is

```@docs
rubisco_limited_rate_glc!
```

Note it here that [`rubisco_limited_rate_glc!`](@ref) only applies to C3
    photosynthesis as the RubisCO-limited rate for C4 plants is
    $V_\text{cmax}$.




## Light-limited Photosynthesis
If leaf internal CO₂ is known, $A_\text{j}$ (gross rate) can be computed using

```@docs
light_limited_rate!
```

If total leaf diffusive conductance to CO₂ is known, $A_\text{j}$ can be
    computed analytically using

```@docs
light_limited_rate_glc!
```

Note it here that [`light_limited_rate_glc!`](@ref) only applies to C3
    photosynthesis as the RubisCO-limited rate for C4 plants is the electron
    transport rate.

Be aware that leaf electron transport rate needs to be calculated before the
    light-limited rate:

```@docs
leaf_ETR!
```




## Product-limited Photosynthesis
If leaf internal CO₂ is known, $A_\text{p}$ (gross rate) can be computed using

```@docs
product_limited_rate!
```

If total leaf diffusive conductance to CO₂ is known, $A_\text{p}$ can be
    computed analytically using

```@docs
product_limited_rate_glc!
```

Note it here that [`product_limited_rate_glc!`](@ref) only applies to C4
    photosynthesis as the RubisCO-limited rate for C4 plants is
    $V_\text{cmax}$/2.




## Photosynthetic Rates
For empirical and optimization stomatal models, iterations are required to get
    the final solution as in StomataModels module. In this case, more
    conveniently computing photosynthetic rates for each leaf is preferable. In
    this case, [`leaf_photo_from_pi!`](@ref) and [`leaf_photo_from_glc!`](@ref)
    are better options:

```@docs
leaf_photo_from_pi!
leaf_photo_from_glc!
```




## Fluorescence
Photosynthesis module also provide ways to compute leaf fluorescence. By
    default, the modules uses fluorescence parameters from Flexas et al. with
    struct [`FluorescenceFlexas`](@ref):

```@docs
AbstractFluoModelParaSet
FluoParaSet
FluorescenceFlexas
```

The function that is used to compute fluorescene is

```@docs
leaf_fluorescence!
```

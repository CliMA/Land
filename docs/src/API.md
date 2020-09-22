# StomtaModels API
```@meta
CurrentModule = StomataModels
```

## Stomatal model schemes
The StomataModels module relies mainly on Photosynthesis and PlantHydraulics
    modules to predict stomatal behavior from plant physiology. This module has
    both empirical and optimal stomatal models. These stomatal models are
    abstractized to an abstract [`AbstractStomatalModel`](@ref), which further
    has subtypes [`EmpiricalStomatalModel`](@ref) and
    [`OptimizationStomatalModel`](@ref).

```@docs
AbstractStomatalModel
EmpiricalStomatalModel
OptimizationStomatalModel
```

Currently, the StomataModels module has four empirical model schemes, and they
    are

```@docs
ESMBallBerry
ESMGentine
ESMLeuning
ESMMedlyn
```

All the empirical models rely on beta functions to make corrections over
    stomatal conductance to account for the stomatal closure with drier soil.
    We have the following prescribed beta function types, and they are:

```@docs
AbstractBetaFunction
```

Some beta functions make correction over the `g1` parameter as in the empitical
    models, and they are:

```@docs
AbstractBetaG
BetaGLinearPleaf
BetaGLinearPsoil
BetaGLinearSWC
```

Some beta functions make correction over the photosynthetic capacity as in the
    Photosynthesis module, and they are:

```@docs
AbstractBetaV
BetaVLinearPleaf
BetaVLinearPsoil
BetaVLinearSWC
```

The beta functions are generalized with

```@docs
β_factor
```

The StomataModels module also contains five optimization model schemes:

```@docs
OSMEller
OSMSperry
OSMWang
OSMWAP
OSMWAPMod
```




## CanopyLayer
The StomataModels module is designed for multi-layer canopies, and each canopy
    has multiple leaves. The stomatal behaviors are modeled per layer basis,
    and the layer may contain any number of leaves starting from 1.
    Photosynthesis-related information is stored in [`CanopyLayer`](@ref)
    struct, but be aware that the leaves have uniform photosynthetic parameters
    and temperature (conductances are different in response to light
    environment).

```@docs
CanopyLayer
```




## Stomatal conductance
For empirical stomatal models, the stomatal conductance is computed as the
    intercept of two functions: an empirical function that describe stomatal
    responses to the physiological and environmental cues and an function that
    follows the diffusion nature of H₂O and CO₂. The abstractized function for
    the empirical correlation is

```@docs
empirical_gsw_from_model
```

For optimization stomatal models, the stomatal conductance is computed as the
    point where the marginal carbon gains equals the marginal carbon risk. The
    marginal carbon gain and risk are generally numerically computed by
    marginally increasing transpiration rate.

This module uses ConstrainedRootSolver module to iterate through the two
    functions to find the solution. The aim is to find the stomatal conductance
    when the [`envir_diff!`](@ref) function equals 0. The [`envir_diff!`](@ref)
    returns the diference between real and model-predicted conductances for
    empirical stomatal models, and the difference between marginal carbon gain
    and risk for optimization stomatal models.

```@docs
envir_diff!
```

In the [`envir_diff!`](@ref) function, leaf photosynthetic rates is modeled
    using [`update_leaf_from_glc!`](@ref), which calculates the gas exchange
    rates from a known total leaf diffusive conductance.

```@docs
update_leaf_from_glc!
```

However, these functions do not force stomatal conductance to stay in its
    ranges. For example, the stomatal conductance solution is set to be zero if
    light is lower than the compensation point. In this case, the
    [`envir_diff!`](@ref) function has to be used along with a control function
    to guarantee realistic stomatal conductance.

```@docs
leaf_gsw_control!
```

To facilitate the use of the StomataModels module, an abstractized function is
    provided for conveniently obtaining stomatal conductance from given
    environmental conditions.

```@docs
leaf_photo_from_envir!
```

To speed up the calculations, leaf physiological parameters are updated only
    if the environmental conditions changes. For example, PAR (photosyntheis
    active radiation) is constant when we iterate [`envir_diff!`](@ref), and
    the electron transport is only updated once. Similar for the cases of
    leaf temperature and soil moisture. This kind of functions used in the
    present module are

```@docs
update_leaf_TP!
update_leaf_AK!
```

I'd like to emphasize it here that the [`leaf_photo_from_envir!`](@ref)
    function only applies to the case of constant leaf temperature because
    leaf energy budget is not calculated, and thus
    [`leaf_photo_from_envir!`](@ref) is only applicable to (1) known leaf
    temperature, and (2) prognostically modeling the non-steady state stomatal
    behaviors. As to the steady state case, leaf energy budget has to be
    considered. For the prognotic stomatal conductance, it is recommended to
    use [`update_leaf_from_gsw!`](@ref) function.

```@docs
update_leaf_from_gsw!
```

Note it here that stomtal conductance is controlled in this function, and thus
    no additional control like [`leaf_gsw_control!`](@ref) is required if
    [`update_leaf_from_gsw!`](@ref) is used.

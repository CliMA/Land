# API
```@meta
CurrentModule = WaterPhysics
```




## Trace molecule
A multiple dispatch approach is used to calculate the temperature and pressure
    dependent physical properties of water and other molecules such as CO₂. The
    trace molecules and mediums are catergorized to `Gas` or `Liquid` subject
    to a general type:

```@docs
AbstractTrace
AbstractTraceGas
AbstractTraceLiquid
```

The supported gasses include

```@docs
TraceGasAir
TraceGasCO₂
TraceGasH₂O
```

The supported liquid includes

```@docs
TraceLiquidH₂O
```




## Capillary pressure
Capillary pressure of liquid water in a pipe is a function of surface tension
    (``\gamma``), pipe raduis (``r``), and contact angle (``\alpha``):

```math
P_{c} = \dfrac{\gamma \cdot \text{cos}(\alpha)}{r}
```

```@docs
capillary_pressure
capillary_pressure(r::FT, T::FT, med::TraceLiquidH₂O{FT}) where {FT<:AbstractFloat}
capillary_pressure(r::FT, T::FT, α::FT, med::TraceLiquidH₂O{FT}) where {FT<:AbstractFloat}
```




## Diffusive coefficient
Diffusion of trace molecules in medium is temperature dependent, to calculate
    this temperature dependency, we provided a function to quick estimate this
    value for different trace molecules

```@docs
diffusive_coefficient
diffusive_coefficient(T::FT, mol::TraceGasCO₂{FT}, med::TraceGasAir{FT}) where {FT<:AbstractFloat}
```

Diffusive coefficient of water vapor is a temperature-dependent function, this
    dependence impact leaf gas exchange via change the maximal stomatal
    conductance to water vapor as the stomatal conductance should not exceed
    its structural limit. To account for this effect, we provided a function to
    calculate the diffusive coefficient relative to 25 Celcius.

```@docs
relative_diffusive_coefficient
relative_diffusive_coefficient(T::FT, mol::AbstractTraceGas{FT}, med::AbstractTraceGas{FT}) where {FT<:AbstractFloat}
```




## Latent heat of evaporation
Water evaporation from liquid phase is a key process to regulate leaf
    temperature, and to best represent this process. We computed the latent
    heat coefficient from water temperature:

```@docs
latent_heat_vapor
latent_heat_vapor(T::FT, med::TraceLiquidH₂O{FT}) where {FT<:AbstractFloat}
```




## Surface tension of air-water interface
When water temperature increases, the surface tension at the air-water
    interface decreases. Surface tension changes impacts the plant water
    transport via two aspects. First, if surface tension is lower, for a
    constant soil water content, the soil matrix potential gets less negative
    because the capillary pressure at the air-water interface decreases. And
    this is beneficial to plants. Second, the air-water interface at the pit
    membrane also has lower capillary pressure when temperature increases,
    meaning that the xylem conduits are less resistant to drought-induced
    air-seeded cavitation. And this is harmful for plants. Though the surface
    tension does not differ much with temperature change within the plant
    physiological active range, we account for this effect in our Land model.

```@docs
surface_tension
surface_tension(T::FT, med::TraceLiquidH₂O{FT}) where {FT<:AbstractFloat}
relative_surface_tension
relative_surface_tension(T::FT, med::TraceLiquidH₂O{FT}) where {FT<:AbstractFloat}
```




## Vapor pressure
When temperature increases, liquid water vapor pressure increases
    exponentially. And this correlation is accounted for using the functions
    below:

```@docs
saturation_vapor_pressure
saturation_vapor_pressure(T::FT, med::TraceLiquidH₂O{FT}) where {FT<:AbstractFloat}
saturation_vapor_pressure(T::FT, Ψ::FT, med::TraceLiquidH₂O{FT}) where {FT<:AbstractFloat}
saturation_vapor_pressure_slope
saturation_vapor_pressure_slope(T::FT, med::TraceLiquidH₂O{FT}) where {FT<:AbstractFloat}
saturation_vapor_pressure_slope(T::FT, Ψ::FT,med::TraceLiquidH₂O{FT}) where {FT<:AbstractFloat}
```

Yet, the saturation vapor pressure is not only a function of temperature, but
    also a function of the air-water interface curvature, known as the Kelvin
    equation. The package provide [`pressure_correction`](@ref) to make the
    correction.

```@docs
pressure_correction
pressure_correction(T::FT, Ψ::FT, med::TraceLiquidH₂O{FT}) where {FT<:AbstractFloat}
```




## Viscosity of liquid water
When temperature increases, liquid water viscosuty decreases, meaning that the
    resistance for water decreases and the pressure drop per flow rate
    decreases. This effect is pretty significant as 1 degree increase of
    temperature results in about 2.4% drop in viscosity, and this is very
    beneficial to plant water transport. Unfortunately, to our knowledge, very
    few models account for this effect when modeling plant hydraulics because
    of the difficulty in modeling the energy budget along the flow path. We
    plan to have this effect accounted for in our CliMA Land model, by
    computing the water tempreature along the flow path and thus the viscosity
    change. The functions to be used are:

```@docs
viscosity
viscosity(T::FT, med::TraceLiquidH₂O{FT}) where {FT<:AbstractFloat}
relative_viscosity
relative_viscosity(T::FT, med::TraceLiquidH₂O{FT}) where {FT<:AbstractFloat}
```

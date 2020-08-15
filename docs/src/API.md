# API
```@meta
CurrentModule = WaterPhysics
```




## Capillary pressure
Capillary pressure of liquid water in a pipe is a function of surface tension
    (``\gamma``), pipe raduis (``r``), and contact angle (``\alpha``):

```math
P_{c} = \dfrac{\gamma \cdot \text{cos}(\alpha)}{r}
```

```@docs
capillary_pressure
```




## Diffusive coefficient of water vapor
Diffusive coefficient of water vapor is a temperature-dependent function, this
    dependence impact leaf gas exchange via change the maximal stomatal
    conductance to water vapor as the stomatal conductance should nto exceed
    its structural limit. To account for this effect, we provided a function to
    calculate the diffusive coefficient relative to 25 Celcius.

```@docs
relative_diffusive_coefficient
```




## Latent heat of evaporation
Water evaporation from liquid phase is a key process to regulate leaf
    temperature, and to best represent this process. We computed the latent
    heat coefficient from water temperature:

```@docs
latent_heat_vapor
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
relative_surface_tension
```




## Vapor pressure
When temperature increases, liquid water vapor pressure increases
    exponentially. And this correlation is accounted for using the functions
    below:

```@docs
saturation_vapor_pressure
saturation_vapor_pressure_slope
```

Yet, the saturation vapor pressure is not only a function of temperature, but
    also a function of the air-water interface curvature, known as the Kelvin
    equation. The package provide [`pressure_correction`](@ref) to make the
    correction.

```@docs
pressure_correction
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
relative_viscosity
```

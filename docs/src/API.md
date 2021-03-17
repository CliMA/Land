# API
```@meta
CurrentModule = WaterPhysics
```




## Trace molecule

```@docs
AbstractTrace
AbstractTraceGas
AbstractTraceLiquid
TraceGasAir
TraceGasCO₂
TraceGasH₂O
TraceGasN₂
TraceGasO₂
TraceLiquidH₂O
```




## Capillary pressure

```@docs
capillary_pressure
capillary_pressure(r::FT, T::FT, med::TraceLiquidH₂O{FT}) where
    {FT<:AbstractFloat}
capillary_pressure(r::FT, T::FT, α::FT, med::TraceLiquidH₂O{FT}) where
    {FT<:AbstractFloat}
```




## Diffusive coefficient

```@docs
diffusive_coefficient
diffusive_coefficient(T::FT, mol::Union{TraceGasCO₂{FT}, TraceGasH₂O{FT},
    TraceGasO₂{FT}}, med::TraceGasAir{FT}) where {FT<:AbstractFloat}
relative_diffusive_coefficient
relative_diffusive_coefficient(T::FT, mol::AbstractTraceGas{FT},
    med::AbstractTraceGas{FT}) where {FT<:AbstractFloat}
relative_diffusive_coefficient(T::FT, mol::Union{TraceGasCO₂{FT},
    TraceGasN₂{FT}}, med::TraceLiquidH₂O{FT}) where {FT<:AbstractFloat}
```




## Latent heat of evaporation

```@docs
latent_heat_vapor
latent_heat_vapor(T::FT, med::TraceLiquidH₂O{FT}) where {FT<:AbstractFloat}
```




## Surface tension of air-water interface

```@docs
surface_tension
surface_tension(T::FT, med::TraceLiquidH₂O{FT}) where {FT<:AbstractFloat}
relative_surface_tension
relative_surface_tension(T::FT, med::TraceLiquidH₂O{FT}) where
    {FT<:AbstractFloat}
```




## Vapor pressure

```@docs
saturation_vapor_pressure
saturation_vapor_pressure(T::FT, med::TraceLiquidH₂O{FT}) where
    {FT<:AbstractFloat}
saturation_vapor_pressure(T::FT, Ψ::FT, med::TraceLiquidH₂O{FT}) where
    {FT<:AbstractFloat}
saturation_vapor_pressure_slope
saturation_vapor_pressure_slope(T::FT, med::TraceLiquidH₂O{FT}) where
    {FT<:AbstractFloat}
saturation_vapor_pressure_slope(T::FT, Ψ::FT,med::TraceLiquidH₂O{FT}) where
    {FT<:AbstractFloat}
pressure_correction
pressure_correction(T::FT, Ψ::FT, med::TraceLiquidH₂O{FT}) where
    {FT<:AbstractFloat}
```




## Viscosity of liquid water

```@docs
viscosity
viscosity(T::FT, med::TraceLiquidH₂O{FT}) where {FT<:AbstractFloat}
relative_viscosity
relative_viscosity(T::FT, med::TraceLiquidH₂O{FT}) where {FT<:AbstractFloat}
```

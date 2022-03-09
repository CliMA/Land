# API
```@meta
CurrentModule = PlantHydraulics
```


## Vulnerability curve
```@docs
relative_hydraulic_conductance
relative_hydraulic_conductance(vc::ComplexVC{FT}, p_25::FT) where {FT<:AbstractFloat}
relative_hydraulic_conductance(vc::LogisticVC{FT}, p_25::FT) where {FT<:AbstractFloat}
relative_hydraulic_conductance(vc::PowerVC{FT}, p_25::FT) where {FT<:AbstractFloat}
relative_hydraulic_conductance(vc::WeibullVC{FT}, p_25::FT) where {FT<:AbstractFloat}
critical_pressure
critical_pressure(vc::ComplexVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat}
critical_pressure(vc::LogisticVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat}
critical_pressure(vc::PowerVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat}
critical_pressure(vc::WeibullVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat}
```

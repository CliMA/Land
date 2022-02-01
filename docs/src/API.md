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
```

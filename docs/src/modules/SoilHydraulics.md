# API
```@meta
CurrentModule = SoilHydraulics
```


## Vulnerability curve
```@docs
BrooksCorey{FT}(vg::VanGenuchten{FT}) where {FT<:AbstractFloat}
soil_θ
soil_ψ_25
relative_hydraulic_conductance
```


## Soil energy+water budget
```@docs
root_sink
soil_budget!
soil_budget!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat}
soil_budget!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, δt::FT) where {FT<:AbstractFloat}
```

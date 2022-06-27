# API
```@meta
CurrentModule = SoilHydraulics
```


## Vulnerability curve
```@docs
BrooksCorey{FT}(vg::VanGenuchten{FT}) where {FT<:AbstractFloat}
soil_θ
soil_θ(bc::BrooksCorey{FT}, ψ_25::FT) where {FT<:AbstractFloat}
soil_θ(vg::VanGenuchten{FT}, ψ_25::FT) where {FT<:AbstractFloat}
soil_ψ_25
soil_ψ_25(bc::BrooksCorey{FT}, θ::FT) where {FT<:AbstractFloat}
soil_ψ_25(vg::VanGenuchten{FT}, θ::FT) where {FT<:AbstractFloat}
relative_hydraulic_conductance
relative_hydraulic_conductance(bc::BrooksCorey{FT}, θ::FT) where {FT<:AbstractFloat}
relative_hydraulic_conductance(bc::BrooksCorey{FT}, ψ::Bool, ψ_25::FT) where {FT<:AbstractFloat}
relative_hydraulic_conductance(vg::VanGenuchten{FT}, θ::FT) where {FT<:AbstractFloat}
relative_hydraulic_conductance(vg::VanGenuchten{FT}, ψ::Bool, ψ_25::FT) where {FT<:AbstractFloat}
```

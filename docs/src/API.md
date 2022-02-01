# API
```@meta
CurrentModule = PlantHydraulics
```


## Vulnerability curve
```@docs
xylem_k_ratio
xylem_k_ratio(vc::LogisticVC{FT}, p_25::FT, vis::FT = FT(1)) where {FT<:AbstractFloat}
xylem_k_ratio(vc::PowerVC{FT}, p_25::FT, vis::FT = FT(1)) where {FT<:AbstractFloat}
xylem_k_ratio(vc::WeibullVC{FT}, p_25::FT, vis::FT = FT(1)) where {FT<:AbstractFloat}
```

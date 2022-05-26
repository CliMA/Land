# API
```@meta
CurrentModule = PlantHydraulics
```


## Vulnerability curve
```@docs
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

## Pressure volume curve
```@docs
xylem_pressure
xylem_pressure(pv::LinearPVCurve{FT}, rvol::FT, T::FT) where {FT<:AbstractFloat}
xylem_pressure(pv::SegmentedPVCurve{FT}, rvol::FT, T::FT) where {FT<:AbstractFloat}
```

## Cavitation legacy
```@docs
clear_legacy!
clear_legacy!(hs::Union{LeafHydraulics{FT}, RootHydraulics{FT}, StemHydraulics{FT}}) where {FT<:AbstractFloat}
clear_legacy!(organ::Union{Leaf{FT}, Root{FT}, Stem{FT}}) where {FT<:AbstractFloat} = clear_legacy!(organ.HS);
clear_legacy!(spac::MonoElementSPAC{FT}) where {FT<:AbstractFloat}
clear_legacy!(spac::MonoGrassSPAC{FT}) where {FT<:AbstractFloat}
clear_legacy!(spac::MonoPalmSPAC{FT}) where {FT<:AbstractFloat}
clear_legacy!(spac::MonoTreeSPAC{FT}) where {FT<:AbstractFloat}
```

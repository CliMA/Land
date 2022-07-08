# API
```@meta
CurrentModule = PlantHydraulics
```


## Vulnerability curve
```@docs
relative_hydraulic_conductance
critical_pressure
```

## Pressure volume curve
```@docs
xylem_pressure
capacitance_buffer
```

## Cavitation legacy
```@docs
clear_legacy!
```

## FLow profile
```@docs
xylem_flow
root_pk
xylem_flow_profile!
xylem_flow_profile!(organ::Union{Leaf{FT}, Leaves2D{FT}, Root{FT}, Stem{FT}}, Δt::FT) where {FT<:AbstractFloat}
xylem_flow_profile!(roots::Vector{Root{FT}}, cache_f::Vector{FT}, cache_k::Vector{FT}, cache_p::Vector{FT}, f_sum::FT, Δt::FT) where {FT<:AbstractFloat}
xylem_flow_profile!(spac::MonoElementSPAC{FT}, Δt::FT) where {FT<:AbstractFloat}
```

## Pressure profile
```@docs
xylem_end_pressure
xylem_pressure_profile!
xylem_pressure_profile!(organ::Union{Leaf{FT}, Leaves2D{FT}, Root{FT}, Stem{FT}}; update::Bool = true) where {FT<:AbstractFloat}
xylem_pressure_profile!(spac::MonoElementSPAC{FT}; update::Bool = true) where {FT<:AbstractFloat}
```

## Critical flow
```@docs
critical_flow
critical_flow(hs::LeafHydraulics{FT}, T::FT, ini::FT = FT(0.5); kr::FT = FT(0.001)) where {FT<:AbstractFloat}
critical_flow(spac::MonoElementSPAC{FT}, ini::FT = FT(0.5); kr::FT = FT(0.001)) where {FT<:AbstractFloat}
```

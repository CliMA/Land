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
capacitance_buffer
capacitance_buffer(pvc::LinearPVCurve{FT}) where {FT<:AbstractFloat}
capacitance_buffer(pvc::SegmentedPVCurve{FT}) where {FT<:AbstractFloat}
```

## Cavitation legacy
```@docs
clear_legacy!
clear_legacy!(hs::Union{LeafHydraulics{FT}, RootHydraulics{FT}, StemHydraulics{FT}}) where {FT<:AbstractFloat}
clear_legacy!(organ::Union{Leaf{FT}, Root{FT}, Stem{FT}}) where {FT<:AbstractFloat}
clear_legacy!(spac::MonoElementSPAC{FT}) where {FT<:AbstractFloat}
clear_legacy!(spac::MonoMLGrassSPAC{FT}) where {FT<:AbstractFloat}
clear_legacy!(spac::MonoMLPalmSPAC{FT}) where {FT<:AbstractFloat}
clear_legacy!(spac::MonoMLTreeSPAC{FT}) where {FT<:AbstractFloat}
```

## FLow profile
```@docs
xylem_flow
xylem_flow(mode::SteadyStateFlow{FT}) where {FT<:AbstractFloat}
xylem_flow(mode::NonSteadyStateFlow{FT}) where {FT<:AbstractFloat}
xylem_flow(hs::Union{LeafHydraulics{FT}, RootHydraulics{FT}, StemHydraulics{FT}}) where {FT<:AbstractFloat}
xylem_flow(organ::Union{Leaf{FT}, Root{FT}, Stem{FT}}) where {FT<:AbstractFloat}
xylem_flow(organs::Vector{Leaf{FT}}) where {FT<:AbstractFloat}
xylem_flow(organs::Vector{Stem{FT}}) where {FT<:AbstractFloat}
root_pk
root_pk(hs::RootHydraulics{FT}, mode::SteadyStateFlow{FT}, T::FT) where {FT<:AbstractFloat}
root_pk(hs::RootHydraulics{FT}, mode::NonSteadyStateFlow{FT}, T::FT) where {FT<:AbstractFloat}
root_pk(hs::RootHydraulics{FT}, T::FT) where {FT<:AbstractFloat}
xylem_flow_profile!
xylem_flow_profile!(mode::SteadyStateFlow, flow::FT) where {FT<:AbstractFloat}
xylem_flow_profile!(mode::NonSteadyStateFlow, f_out::FT) where {FT<:AbstractFloat}
xylem_flow_profile!(hs::Union{LeafHydraulics{FT}, RootHydraulics{FT}, StemHydraulics{FT}}, mode::SteadyStateFlow{FT}, T::FT, Δt::FT) where {FT<:AbstractFloat}
xylem_flow_profile!(hs::LeafHydraulics{FT}, mode::NonSteadyStateFlow{FT}, T::FT, Δt::FT) where {FT<:AbstractFloat}
xylem_flow_profile!(hs::Union{RootHydraulics{FT}, StemHydraulics{FT}}, mode::NonSteadyStateFlow{FT}, T::FT, Δt::FT) where {FT<:AbstractFloat}
xylem_flow_profile!(hs::Union{LeafHydraulics{FT}, RootHydraulics{FT}, StemHydraulics{FT}}, T::FT, Δt::FT) where {FT<:AbstractFloat}
xylem_flow_profile!(organ::Union{Leaf{FT}, Root{FT}, Stem{FT}}, Δt::FT) where {FT<:AbstractFloat}
xylem_flow_profile!(roots::Vector{Root{FT}}, cache_f::Vector{FT}, cache_k::Vector{FT}, cache_p::Vector{FT}, f_sum::FT, Δt::FT) where {FT<:AbstractFloat}
xylem_flow_profile!(spac::MonoElementSPAC{FT}, Δt::FT) where {FT<:AbstractFloat}
xylem_flow_profile!(spac::MonoMLGrassSPAC{FT}, Δt::FT) where {FT<:AbstractFloat}
xylem_flow_profile!(spac::MonoMLPalmSPAC{FT}, Δt::FT) where {FT<:AbstractFloat}
xylem_flow_profile!(spac::MonoMLTreeSPAC{FT}, Δt::FT) where {FT<:AbstractFloat}
```

## Pressure profile
```@docs
xylem_end_pressure
xylem_end_pressure(hs::LeafHydraulics{FT}, flow::FT, T::FT) where {FT<:AbstractFloat}
xylem_end_pressure(hs::RootHydraulics{FT}, flow::FT, T::FT) where {FT<:AbstractFloat}
xylem_end_pressure(hs::StemHydraulics{FT}, flow::FT, T::FT) where {FT<:AbstractFloat}
xylem_end_pressure(organ::Union{Leaf{FT}, Root{FT}, Stem{FT}}, flow::FT) where {FT<:AbstractFloat}
xylem_end_pressure(spac::MonoElementSPAC{FT}, flow::FT) where {FT<:AbstractFloat}
xylem_end_pressure(spac::MonoElementSPAC{FT}, f_sl::FT, f_sh::FT, r_sl::FT) where {FT<:AbstractFloat}
xylem_pressure_profile!
xylem_pressure_profile!(hs::LeafHydraulics{FT}, mode::SteadyStateFlow{FT}, T::FT; update::Bool = true) where {FT<:AbstractFloat}
xylem_pressure_profile!(hs::LeafHydraulics{FT}, mode::NonSteadyStateFlow{FT}, T::FT; update::Bool = true) where {FT<:AbstractFloat}
xylem_pressure_profile!(hs::RootHydraulics{FT}, mode::SteadyStateFlow{FT}, T::FT; update::Bool = true) where {FT<:AbstractFloat}
xylem_pressure_profile!(hs::RootHydraulics{FT}, mode::NonSteadyStateFlow{FT}, T::FT; update::Bool = true) where {FT<:AbstractFloat}
xylem_pressure_profile!(hs::StemHydraulics{FT}, mode::SteadyStateFlow{FT}, T::FT; update::Bool = true) where {FT<:AbstractFloat}
xylem_pressure_profile!(hs::StemHydraulics{FT}, mode::NonSteadyStateFlow{FT}, T::FT; update::Bool = true) where {FT<:AbstractFloat}
xylem_pressure_profile!(organ::Union{Leaf{FT}, Root{FT}, Stem{FT}}; update::Bool = true) where {FT<:AbstractFloat}
xylem_pressure_profile!(spac::MonoElementSPAC{FT}; update::Bool = true) where {FT<:AbstractFloat}
xylem_pressure_profile!(spac::MonoMLGrassSPAC{FT}; update::Bool = true) where {FT<:AbstractFloat}
xylem_pressure_profile!(spac::MonoMLPalmSPAC{FT}; update::Bool = true) where {FT<:AbstractFloat}
xylem_pressure_profile!(spac::MonoMLTreeSPAC{FT}; update::Bool = true) where {FT<:AbstractFloat}
```

## Critical flow
```@docs
critical_flow
critical_flow(hs::LeafHydraulics{FT}, T::FT, ini::FT = FT(0.5); kr::FT = FT(0.001)) where {FT<:AbstractFloat}
critical_flow(spac::MonoElementSPAC{FT}, ini::FT = FT(0.5); kr::FT = FT(0.001)) where {FT<:AbstractFloat}
```

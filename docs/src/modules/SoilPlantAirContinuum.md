# SoilPlantAirContinuum
```@meta
CurrentModule = SoilPlantAirContinuum
```

## Core function
```@docs
initialize!
soil_plant_air_continuum!
soil_plant_air_continuum!(spac::Union{MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC{FT}}, δt::FT; update::Bool = false, θ_on::Bool = true) where {FT<:AbstractFloat}
adjusted_time
time_stepper!
```

## Update SPAC
```@docs
update!
update!(air::AirLayer{FT}; p_CO₂::Union{Number,Nothing} = nothing, p_H₂O::Union{Number,Nothing} = nothing, rh::Union{Number,Nothing} = nothing, t::Union{Number,Nothing} = nothing, vpd::Union{Number,Nothing} = nothing, wind::Union{Number,Nothing} = nothing) where {FT<:AbstractFloat}
```

## Measures
```@docs
CNPP
GPP
PPAR
T_VEG
```

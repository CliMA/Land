module WaterPhysics

using Parameters

using ..LandParameters

@unpack CP_I,
        CP_L,
        CP_V,
        K_25,
        LH_S0,
        LH_V0,
        PRESS_TRIPLE,
        R_V,
        T_TRIPLE = LandParameters

export latent_heat_vapor,
       relative_diffusive_coefficient,
       relative_surface_tension,
       relative_viscosity,
       saturation_vapor_pressure,
       saturation_vapor_pressure_slope,
       surface_tension,
       viscosity




###############################################################################
#
# Latent heat
# Adapted from the ClimateMachine/Common/ThermoDynamics/relations.jl
# These functions passed the FT test
# These functions are documented in the Utils page
#
###############################################################################
"""
    latent_heat_vapor(T::FT)
    latent_heat_vapor(T::Array)

The specific latent heat of vaporization, given
 - `T` (Array of) temperature
"""
function latent_heat_vapor(T::FT) where {FT <: Real}
    _cp_i ::FT = CP_I
    _cp_l ::FT = CP_L
    _cp_v ::FT = CP_V
    _LH_s0::FT = LH_S0
    _LH_v0::FT = LH_V0
    _T_0  ::FT = T_TRIPLE
    _Δcp_s::FT = _cp_v - _cp_i
    _Δcp_v::FT = _cp_v - _cp_l
    if T >= _T_0
        return _LH_v0 + _Δcp_v * (T - _T_0)
    else
        return _LH_s0 + _Δcp_s * (T - _T_0)
    end
end

function latent_heat_vapor(T::Array)
    FT         = eltype(T)
    _cp_i ::FT = CP_I
    _cp_l ::FT = CP_L
    _cp_v ::FT = CP_V
    _LH_s0::FT = LH_S0
    _LH_v0::FT = LH_V0
    _T_0  ::FT = T_TRIPLE
    _Δcp_s::FT = _cp_v - _cp_i
    _Δcp_v::FT = _cp_v - _cp_l

    _LH = T .* 0
    mask = (T .>= _T_0)
    _LH[  mask] .= _LH_v0 .+ _Δcp_v .* (T[  mask] .- _T_0)
    _LH[.!mask] .= _LH_v0 .+ _Δcp_v .* (T[.!mask] .- _T_0)
    return _LH
end








###############################################################################
#
# Saturated vapor pressure
# Adapted from the ClimateMachine/Common/ThermoDynamics/relations.jl
# TODO Merge into ThermoDynamics?
# These functions passed the FT test
# These functions are documented in the Utils page
#
###############################################################################
"""
    saturation_vapor_pressure(T::FT)
    saturation_vapor_pressure(T::Array)

Return the saturation vapor pressure over a plane liquid surface given
 - `T` (Array of) temperature
"""
function saturation_vapor_pressure(T::FT) where {FT}
    _LH_v0       ::FT = LH_V0
    _cp_v        ::FT = CP_V
    _cp_l        ::FT = CP_L
    _press_triple::FT = PRESS_TRIPLE
    _R_v         ::FT = R_V
    _T_triple    ::FT = T_TRIPLE
    _T_0         ::FT = T_TRIPLE
    _Δcp         ::FT = _cp_v - _cp_l

    return _press_triple *
           (T / _T_triple)^(_Δcp / _R_v) *
           exp((_LH_v0 - _Δcp * _T_0) / _R_v * (1 / _T_triple - 1 / T))
end

function saturation_vapor_pressure(T::Array)
    FT = eltype(T)

    _LH_v0       ::FT = LH_V0
    _cp_v        ::FT = CP_V
    _cp_l        ::FT = CP_L
    _press_triple::FT = PRESS_TRIPLE
    _R_v         ::FT = R_V
    _T_triple    ::FT = T_TRIPLE
    _T_0         ::FT = T_TRIPLE
    _Δcp         ::FT = _cp_v - _cp_l

    return _press_triple .*
           (T ./ _T_triple) .^ (_Δcp / _R_v) .*
           exp.((_LH_v0 - _Δcp * _T_0) ./ _R_v .* (1 / _T_triple .- 1 ./ T))
end




"""
    saturation_vapor_pressure_slope(T::FT)
    saturation_vapor_pressure_slope(T::Array)

Return the 1st order derivative of saturation vapor pressure over a plane liquid surface given
 - `T` (Array of) temperature

Compute the the 1st order derivative of saturation vapor pressure over a plane surface by integration of the Clausius-Clapeyron relation.

The re-arranged Clausius-Clapeyron relation: `dp_v_sat/dT = [LH_0 + Δcp * (T-T_triple)]/(R_v*T^2) * p_v_sat`.
"""
function saturation_vapor_pressure_slope(T::FT) where {FT}
    _LH_v0::FT = LH_V0
    _cp_v ::FT = CP_V
    _cp_l ::FT = CP_L
    _R_v  ::FT = R_V
    _T_0  ::FT = T_TRIPLE
    _Δcp  ::FT = _cp_v - _cp_l

    return (_LH_v0 + _Δcp * (T - _T_0)) / (_R_v*T^2)
end

function saturation_vapor_pressure_slope(T::Array)
    FT         = eltype(T)
    _LH_v0::FT = LH_V0
    _cp_v ::FT = CP_V
    _cp_l ::FT = CP_L
    _R_v  ::FT = R_V
    _T_0  ::FT = T_TRIPLE
    _Δcp  ::FT = _cp_v - _cp_l

    return (_LH_v0 .+ _Δcp .* (T .- _T_0)) ./ (_R_v .* T.^2)
end








###############################################################################
#
# Surface tension functions
# These functions passed the FT test
# These functions are documented in the Utils page
#
###############################################################################
const ST_k      = 0.2358
const ST_T_crit = 647.096
const ST_exp    = 1.256
const ST_corr   = 0.625
const ST_ref    = 0.07197220523

"""
    surface_tension(T::FT)
    surface_tension(T::Array)

Surface tension `[N m⁻¹]` of water against air, given
- `T` (Array of) temperature

The equation used is `γ = 0.2358 * (1 - T/T_c)^1.256 * [1 - 0.625 * (1 - T/T_c)]`
See http://www.iapws.org/relguide/Surf-H2O.html
"""
function surface_tension(T::FT) where {FT}
    _ST_corr    ::FT = ST_corr
    _ST_exp     ::FT = ST_exp
    _ST_k       ::FT = ST_k
    _ST_T_crit  ::FT = ST_T_crit
    _ST_T_r_diff::FT = 1 - T / _ST_T_crit
    return _ST_k * _ST_T_r_diff^_ST_exp * (1 - _ST_corr*_ST_T_r_diff)
end

function surface_tension(T::Array)
    FT               = eltype(T)
    _ST_corr    ::FT = ST_corr
    _ST_exp     ::FT = ST_exp
    _ST_k       ::FT = ST_k
    _ST_T_crit  ::FT = ST_T_crit
    _ST_T_r_diff     = 1 .- T ./ _ST_T_crit
    return _ST_k .* _ST_T_r_diff .^ _ST_exp .* (1 .- _ST_corr .* _ST_T_r_diff)
end




"""
    relative_surface_tension(T::FT)
    relative_surface_tension(T::Array)

Relative surface tension of water against air relative to 298.15 K, given
- `T` (Array of) temperature

The equation used is `γ = 0.2358 * (1 - T/T_c)^1.256 * [1 - 0.625 * (1 - T/T_c)]`.
See http://www.iapws.org/relguide/Surf-H2O.html

A reference temperature at 298.15 K is used because the hydraulic conductances and vulnerability curves of plants are described at 298.15 K.
"""
function relative_surface_tension(T::FT) where {FT}
    _ST_ref::FT = ST_ref
    return surface_tension(T) / _ST_ref
end

function relative_surface_tension(T::Array)
    FT = eltype(T)
    _ST_ref::FT = ST_ref
    return surface_tension(T) ./ _ST_ref
end








###############################################################################
#
# Viscosity functions
# These functions passed the FT test
# These functions are documented in the Utils page
#
###############################################################################
const VIS_A =  1.856e-14    # Pa s | Used for viscosity
const VIS_B =  4209.0       # K    | Used for viscosity
const VIS_C =  0.04527      # K⁻¹  | Used for viscosity
const VIS_D = -3.376e-5     # K⁻²  | Used for viscosity

"""
    viscosity(T::FT)
    viscosity(T::Array)

Viscosity of water `[Pa s]`, given
- `T` (Array of) temperature

Equation used is `υ = A * exp( B/T + C*T + D*T^2 )`, and the fitting parameters are from Reid, Prausnitz, & Poling (1987), valid through 273-643 K
```
A = 1.856E-14 # Pa s
B = 4209      # K
C = 0.04527   # K⁻¹
D = -3.376E-5 # K⁻²
```
"""
function viscosity(T::FT) where {FT}
    return FT(VIS_A) * exp(FT(VIS_B)/T + FT(VIS_C)*T + FT(VIS_D)*T^2)
end

function viscosity(T::Array)
    FT = eltype(T)
    return FT(VIS_A) .* exp.( FT(VIS_B) ./ T .+ FT(VIS_C) .* T .+ FT(VIS_D) .* T.^2 )
end




"""
    relative_viscosity(T::FT)
    relative_viscosity(T::Array)

Viscosity relative to 25 degree C (298.15 K), given
- `T` Water temperature

Equation used is `υ/υ25 = exp( B/T + C*T + D*T^2 - B/tem25 - C*tem25 - D*tem25^2 )`, and the fitting parameters are from Reid, Prausnitz, & Poling (1987), valid through 273-643 K
```
B = 4209      # K
C = 0.04527   # K⁻¹
D = -3.376E-5 # K⁻²
```
"""
function relative_viscosity(T::FT) where {FT}
    return exp( FT(VIS_B) * ( 1/T - 1/FT(K_25)) + FT(VIS_C) * (T - FT(K_25)) + FT(VIS_D) * (T^2 - FT(K_25)^2) )
end

function relative_viscosity(T::Array)
    FT = eltype(T)
    return exp.( FT(VIS_B) .* ( 1 ./ T .- 1/FT(K_25)) .+ FT(VIS_C) .* (T .- FT(K_25)) .+ FT(VIS_D) .* (T.^2 .- FT(K_25)^2) )
end








###############################################################################
#
# Functions in development
# Test and document the functions before merging them to this file
#
###############################################################################
include("WaterPhysics_in_development.jl")




end # module

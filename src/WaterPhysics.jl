module WaterPhysics

using CLIMAParameters
using CLIMAParameters.Planet: LH_v0, R_v, T_freeze, T_triple, cp_l, cp_v,
            molmass_water, press_triple, ρ_cloud_liq
using DocStringExtensions: METHODLIST, TYPEDEF, TYPEDFIELDS
using Parameters: @unpack




# Define a local struct inherited from AbstractEarthParameterSet
struct EarthParameterSet <: AbstractEarthParameterSet end
const EARTH      = EarthParameterSet();
CP_L(FT)         = FT( cp_l(EARTH) );
CP_V(FT)         = FT( cp_v(EARTH) );
GAS_R(FT)        = FT( gas_constant() );
K_25(FT)         = FT( T_freeze(EARTH) + 25 );
LH_V0(FT)        = FT( LH_v0(EARTH) );
MOLMASS_H₂O(FT)  = FT( molmass_water(EARTH) );
PRESS_TRIPLE(FT) = FT( press_triple(EARTH) );
R_V(FT)          = FT( R_v(EARTH) );
T_TRIPLE(FT)     = FT( T_triple(EARTH) );
ρ_H₂O(FT)        = FT( ρ_cloud_liq(EARTH) );

MOLVOL_H₂O(FT)   = MOLMASS_H₂O(FT) / ρ_H₂O(FT);




# export public types
export TraceGasAir,
       TraceGasCO₂,
       TraceGasH₂O,
       TraceLiquidH₂O

# export public functions
export capillary_pressure,
       diffusive_coefficient,
       latent_heat_vapor,
       pressure_correction,
       relative_diffusive_coefficient,
       relative_surface_tension,
       relative_viscosity,
       saturation_vapor_pressure,
       saturation_vapor_pressure_slope,
       surface_tension,
       viscosity




###############################################################################
#
# Trace gas or liquid
#
###############################################################################
"""
Identity label for trace gas or liquid. Hierarchy of `AbstractTrace`:
- [`AbstractTraceGas`](@ref)
- [`AbstractTraceLiquid`](@ref)

$(TYPEDEF)

"""
abstract type AbstractTrace{FT<:AbstractFloat} end




"""
Identity label for trace gas. Hierarchy of `AbstractTraceGas`:
- [`TraceGasAir`](@ref)
- [`TraceGasCO₂`](@ref)
- [`TraceGasH₂O`](@ref)

$(TYPEDEF)

"""
abstract type AbstractTraceGas{FT<:AbstractFloat} <: AbstractTrace{FT} end




"""
Identity label for trace liquid. Hierarchy of `AbstractTraceLiquid`:
- [`TraceLiquidH₂O`](@ref)

$(TYPEDEF)

"""
abstract type AbstractTraceLiquid{FT<:AbstractFloat} <: AbstractTrace{FT} end




"""
Identity label for air.

$(TYPEDEF)

"""
struct TraceGasAir{FT<:AbstractFloat} <: AbstractTraceGas{FT} end




"""
Identity label for gas phase CO₂.

$(TYPEDEF)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct TraceGasCO₂{FT<:AbstractFloat} <: AbstractTraceGas{FT}
    "Diffusive coefficient in air"
    d_air::FT = 1.6e-5
end




"""
Identity label for gas phase H₂O.

$(TYPEDEF)

"""
struct TraceGasH₂O{FT<:AbstractFloat} <: AbstractTraceGas{FT} end




"""
Identity label for liquid phase H₂O.

$(TYPEDEF)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct TraceLiquidH₂O{FT<:AbstractFloat} <: AbstractTraceLiquid{FT}
    # Surface tension related
    # γ = γ_k * (1-T/γ_T_c)^γ_exp * [1 - γ_cor*(1-T/γ_T_c)]
    "Surface tension coefficient correction"
    γ_cor::FT = 0.625
    "Surface tension coefficient exponent"
    γ_exp::FT = 1.256
    "Surface tension coefficient k in `[N m⁻¹]`"
    γ_k  ::FT = 0.2358
    "Surface tension at 298.15 K in `[N m⁻¹]`"
    γ_ref::FT = 0.07197220523
    "Surface tension coefficient T_crit in `[K]`"
    γ_T_c::FT = 647.096

    # Viscosity related
    # υ = A ⋅ exp( A/T + C⋅T + D⋅T^2 )
    "Viscosity coefficient A in `[Pa s]`"
    υ_A::FT = 1.856e-14
    "Viscosity coefficient B in `[K]`"
    υ_B::FT = 4209.0
    "Viscosity coefficient C in `[K⁻¹]`"
    υ_C::FT = 0.04527
    "Viscosity coefficient D in `[K⁻²]`"
    υ_D::FT = -3.376e-5
end








###############################################################################
#
# Capillary pressure
#
###############################################################################
"""
Capillary pressure of trace liquid.

$(METHODLIST)

"""
function capillary_pressure end




"""
    capillary_pressure(
                r::FT,
                T::FT,
                med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
    ) where {FT<:AbstractFloat}

Capillary pressure of trace liquid in `[Pa]`, given
- `r` Curvature radius in `[m]`
- `T` Trace liquid temperature in `[K]`
- `med` Medium. Optional. Default is liquid water
"""
capillary_pressure(
            r::FT,
            T::FT,
            med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
(
     return 2 * surface_tension(T, med) / r
)




"""
    capillary_pressure(
                r::FT,
                T::FT,
                α::FT,
                med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
    ) where {FT<:AbstractFloat}

Capillary pressure of trace liquid in `[Pa]`, given
- `r` Curvature radius in `[m]`
- `T` Trace liquid temperature in `[K]`
- `α` Contact angle in `[°]`
- `med` Medium. Optional. Default is liquid water
"""
capillary_pressure(
            r::FT,
            T::FT,
            α::FT,
            med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
(
     return cosd(α) * capillary_pressure(r, T, med)
)








###############################################################################
#
# Diffusive coefficient
#
###############################################################################
"""
Diffusive coefficient of trace molecule in medium

$(METHODLIST)

"""
function diffusive_coefficient end




"""
    diffusive_coefficient(
                T::FT,
                mol::TraceGasCO₂{FT},
                med::TraceGasAir{FT}
    ) where {FT<:AbstractFloat}

Diffusive coefficient of trace molecule in `[m² s⁻¹]`, given
- `T` Trace medium temperature in `[K]`
- `mol` Trace molecule
- `med` Diffusion medium
"""
diffusive_coefficient(
            T::FT,
            mol::TraceGasCO₂{FT},
            med::TraceGasAir{FT}
) where {FT<:AbstractFloat} =
(
    return mol.d_air * relative_diffusive_coefficient(T, mol, med)
)




"""
Relative diffusive coefficient of trace gas in medium

$(METHODLIST)

"""
function relative_diffusive_coefficient end




"""
    relative_diffusive_coefficient(
                T::FT,
                mol::AbstractTraceGas{FT} = TraceGasH₂O{FT}(),
                med::AbstractTraceGas{FT} = TraceGasAir{FT}()
    ) where {FT<:AbstractFloat}

Relative diffusive coefficient of trace gas in medium, given
- `T` Water vapor temperature in `[K]`
- `mol` Trace molecule. Optional, default is water vapor
- `med` Medium. Optional, default is air
"""
relative_diffusive_coefficient(
            T::FT,
            mol::AbstractTraceGas{FT} = TraceGasH₂O{FT}(),
            med::AbstractTraceGas{FT} = TraceGasAir{FT}()
) where {FT<:AbstractFloat} =
(
    return (T / K_25(FT)) ^ FT(1.8)
)








###############################################################################
#
# Latent heat of evaporation
# Adapted from the ClimateMachine/Common/ThermoDynamics/relations.jl
#
###############################################################################
"""
Latent heat of vaporization.

$(METHODLIST)

"""
function latent_heat_vapor end




"""
    latent_heat_vapor(
                T::FT,
                med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
    ) where {FT<:AbstractFloat}

Latent heat of vaporization in `[J kg⁻¹]`, given
- `T` Medium temperature in `[K]`
- `med` Medium. Optional. Default is liquid water
"""
latent_heat_vapor(
            T::FT,
            med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
(
    _LH_v0::FT = LH_V0(FT);
    _T_0  ::FT = T_TRIPLE(FT);
    _Δcp_v::FT = CP_V(FT) - CP_L(FT);

    return _LH_v0 + _Δcp_v * (T - _T_0)
)








###############################################################################
#
# Saturated vapor pressure
# Adapted from the ClimateMachine/Common/ThermoDynamics/relations.jl
#
###############################################################################
"""
Kelvin correction factor for saturation vapor pressure.

$(METHODLIST)

"""
function pressure_correction end




"""
    pressure_correction(
                T::FT,
                Ψ::FT,
                med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
    ) where {FT<:AbstractFloat}

Kelvin correction factor for saturation vapor pressure, given
- `T` Liquid water temperature in `[K]`
- `Ψ` Liquid water pressure in `[Pa]`, positive/negative for convex/concave
    interface
- `med` Medium. Optional. Default is liquid water

The Kelvin equation is
```math
\\log \\left( \\dfrac{P_{sat}}{P_{sat}^{*}} \\right) =
\\dfrac{Ψ ⋅ V_{m}}{R ⋅ T}
```
"""
pressure_correction(
            T::FT,
            Ψ::FT,
            med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
(
    return exp((Ψ * MOLVOL_H₂O(FT)) / (GAS_R(FT) * T))
)




"""
Saturation vapor pressure.

$(METHODLIST)

"""
function saturation_vapor_pressure end




"""
    saturation_vapor_pressure(
                T::FT,
                med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
    ) where {FT<:AbstractFloat}

Saturation vapor pressure in `[Pa]`, given
- `T` Liquid water temperature in `[K]`
- `med` Medium. Optional. Default is liquid water
"""
saturation_vapor_pressure(
            T::FT,
            med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
(
    _LH_v0       ::FT = LH_V0(FT);
    _cp_v        ::FT = CP_V(FT);
    _cp_l        ::FT = CP_L(FT);
    _press_triple::FT = PRESS_TRIPLE(FT);
    _R_v         ::FT = R_V(FT);
    _T_triple    ::FT = T_TRIPLE(FT);
    _T_0         ::FT = T_TRIPLE(FT);
    _Δcp         ::FT = _cp_v - _cp_l;

    return _press_triple *
           (T / _T_triple)^(_Δcp / _R_v) *
           exp((_LH_v0 - _Δcp * _T_0) / _R_v * (1 / _T_triple - 1 / T))
)




"""
    saturation_vapor_pressure(
                T::FT,
                Ψ::FT,
                med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
    ) where {FT<:AbstractFloat}

Saturation vapor pressure in `[Pa]`, given
- `T` Liquid water temperature in `[K]`
- `Ψ` Liquid water pressure in `[Pa]`, positive/negative for convex/concave
    interface; if `Ψ` is given, [`pressure_correction`](@ref) is made
- `med` Medium. Optional. Default is liquid water
"""
saturation_vapor_pressure(
            T::FT,
            Ψ::FT,
            med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
(
    return saturation_vapor_pressure(T, med) * pressure_correction(T, Ψ, med)
)




"""
First order derivative of saturation vapor pressure.

$(METHODLIST)

"""
function saturation_vapor_pressure_slope end




"""
    saturation_vapor_pressure_slope(
                T::FT,
                med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
    ) where {FT<:AbstractFloat}

First order derivative of saturation vapor pressure in `[Pa K⁻¹]`, given
- `T` Liquid water temperature in `[K]`
- `med` Medium. Optional. Default is liquid water

Compute the the 1st order derivative of saturation vapor pressure over a plane
    surface by integration of the Clausius-Clapeyron relation.

The re-arranged Clausius-Clapeyron relation
```math
\\frac{∂P_{sat}^{*}}{∂T} = P_{sat} ⋅ \\dfrac{ LH_0 + Δc_p ⋅ (T - T_{triple})}
                                        { R_v ⋅ T^2 }
```
"""
saturation_vapor_pressure_slope(
            T::FT,
            med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
(
    _LH_v0::FT = LH_V0(FT);
    _cp_v ::FT = CP_V(FT);
    _cp_l ::FT = CP_L(FT);
    _R_v  ::FT = R_V(FT);
    _T_0  ::FT = T_TRIPLE(FT);
    _Δcp  ::FT = _cp_v - _cp_l;

    return (_LH_v0 + _Δcp * (T - _T_0)) / (_R_v*T^2)
)




"""
    saturation_vapor_pressure_slope(
                T::FT,
                Ψ::FT,
                med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
    ) where {FT<:AbstractFloat}

First order derivative of saturation vapor pressure in `[Pa K⁻¹]`, given
- `T` Liquid water temperature in `[K]`
- `Ψ` Liquid water pressure in `[Pa]`, positive/negative for convex/concave
    interface; if `Ψ` is given, [`pressure_correction`](@ref) is made
- `med` Medium. Optional. Default is liquid water
"""
saturation_vapor_pressure_slope(
            T::FT,
            Ψ::FT,
            med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
(
    return saturation_vapor_pressure_slope(T, med) *
           pressure_correction(T, Ψ, med)
)








###############################################################################
#
# Surface tension
#
###############################################################################
"""
Surface tension.

$(METHODLIST)

"""
function surface_tension end




"""
    surface_tension(
                T::FT,
                med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
    ) where {FT<:AbstractFloat}

Surface tension of trace liquid in `[N m⁻¹]`, given
- `T` Liquid water temperature in `[K]`
- `med` Medium. Optional. Default is liquid water

The equation used is
```math
γ = 0.2358 ⋅
    \\left( 1 - \\dfrac{T}{T_c} \\right)^{1.256} ⋅
    \\left[ 1 - 0.625 ⋅ \\left( 1 - \\dfrac{T}{T_c} \\right) \\right]
```
See http://www.iapws.org/relguide/Surf-H2O.html
"""
surface_tension(
            T::FT,
            med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
(
    @unpack γ_cor, γ_exp, γ_k, γ_T_c = med;
    _γ_T_r_diff = 1 - T / γ_T_c;

    return γ_k * _γ_T_r_diff^γ_exp * (1 - γ_cor*_γ_T_r_diff)
)




"""
Relative surface tension of trace liquid relative to 298.15 K.

$(METHODLIST)

"""
function relative_surface_tension end




"""
    relative_surface_tension(
                T::FT,
                med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
    ) where {FT<:AbstractFloat}

Relative surface tension of trace liquid relative to 298.15 K, given
- `T` Liquid water temperature in `[K]`
- `med` Medium. Optional. Default is liquid water
"""
relative_surface_tension(
            T::FT,
            med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
(
    return surface_tension(T, med) / med.γ_ref
)








###############################################################################
#
# Viscosity
#
###############################################################################
"""
Viscosity.

$(METHODLIST)

"""
function viscosity end




"""
    viscosity(T::FT,
              med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
    ) where {FT<:AbstractFloat}

Viscosity in `[Pa s]`, given
- `T` Liquid water temperature in `[K]`
- `med` Medium. Optional. Default is liquid water

Equation used is
```math
υ = A ⋅ \\exp \\left( \\dfrac{B}{T} + C⋅T + D⋅T^2 \\right)
```
and the fitting parameters are from Reid, Prausnitz, & Poling (1987), valid
through 273-643 K
```
A = 1.856E-14 # Pa s
B = 4209      # K
C = 0.04527   # K⁻¹
D = -3.376E-5 # K⁻²
```
"""
viscosity(T::FT,
          med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
(
    @unpack υ_A, υ_B, υ_C, υ_D = med;

    return υ_A * exp( υ_B/T + υ_C*T + υ_D*T^2 )
)




"""
Viscosity relative to 298.15 K

$(METHODLIST)

"""
function relative_viscosity end




"""
    relative_viscosity(
                T::FT,
                med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
    ) where {FT<:AbstractFloat}

Viscosity relative to 298.15 K, given
- `T` Liquid water temperature in `[K]`
- `liquid` Optional. Default is liquid water

Equation used is
```math
υ = A ⋅ \\exp \\left( \\dfrac{B}{T} + C⋅T + D⋅T^2 \\right)
```
and the fitting parameters are from Reid, Prausnitz, & Poling (1987), valid
through 273-643 K
```
B = 4209      # K
C = 0.04527   # K⁻¹
D = -3.376E-5 # K⁻²
```
"""
relative_viscosity(
            T::FT,
            med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
(
    @unpack υ_B, υ_C, υ_D = med;
    _K = K_25(FT);

    return exp( υ_B * ( 1/T - 1/_K) + υ_C * (T - _K) + υ_D * (T^2 - _K^2) )
)



end # module

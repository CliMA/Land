module WaterPhysics

using CLIMAParameters
using CLIMAParameters.Planet: LH_v0, R_v, T_freeze, T_triple, cp_l, cp_v,
            molmass_water, press_triple, ρ_cloud_liq
using DocStringExtensions: METHODLIST, TYPEDEF




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
export AbstractTraceGas,
       AbstractTraceLiquid,
       AbstractTraceMolecule,
       TraceGasAir,
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
# Abstracted trace gas
#
###############################################################################
"""
Identity label for trace gas or liquid. Hierarchy of `AbstractTraceMolecule`:
- [`AbstractTraceGas`](@ref)
- [`AbstractTraceLiquid`](@ref)

$(TYPEDEF)

"""
abstract type AbstractTraceMolecule end




"""
Identity label for trace gas. Hierarchy of `AbstractTraceGas`:
- [`TraceGasAir`](@ref)
- [`TraceGasCO₂`](@ref)
- [`TraceGasH₂O`](@ref)

$(TYPEDEF)

"""
abstract type AbstractTraceGas <: AbstractTraceMolecule end




"""
Identity label for trace liquid. Hierarchy of `AbstractTraceLiquid`:
- [`TraceLiquidH₂O`](@ref)

$(TYPEDEF)

"""
abstract type AbstractTraceLiquid <: AbstractTraceMolecule end




"""
Identity label for air.

$(TYPEDEF)

"""
struct TraceGasAir <: AbstractTraceGas end




"""
Identity label for gas phase CO₂.

$(TYPEDEF)

"""
struct TraceGasCO₂ <: AbstractTraceGas end




"""
Identity label for gas phase H₂O.

$(TYPEDEF)

"""
struct TraceGasH₂O <: AbstractTraceGas end




"""
Identity label for liquid phase H₂O.

$(TYPEDEF)

"""
struct TraceLiquidH₂O <: AbstractTraceLiquid end








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
                med::TraceLiquidH₂O = TraceLiquidH₂O()
    ) where {FT<:AbstractFloat}

Capillary pressure of trace liquid in `[Pa]`, given
- `r` Curvature radius in `[m]`
- `T` Trace liquid temperature in `[K]`
- `med` Medium. Optional. Default is liquid water
"""
capillary_pressure(
            r::FT,
            T::FT,
            med::TraceLiquidH₂O = TraceLiquidH₂O()
) where {FT<:AbstractFloat} =
(
     return 2 * surface_tension(T, med) / r
)




"""
    capillary_pressure(
                r::FT,
                T::FT,
                α::FT,
                med::TraceLiquidH₂O = TraceLiquidH₂O()
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
            med::TraceLiquidH₂O = TraceLiquidH₂O()
) where {FT<:AbstractFloat} =
(
     return cosd(α) * capillary_pressure(r, T, med)
)








###############################################################################
#
# Diffusive coefficient vs T
#
###############################################################################
const D_CO2_IN_AIR = 1.6e-5;
"""
Diffusive coefficient of trace molecule in medium

$(METHODLIST)

"""
function diffusive_coefficient end




"""
    diffusive_coefficient(
                T::FT,
                mol::TraceGasCO₂,
                med::TraceGasAir
    ) where {FT<:AbstractFloat}

Diffusive coefficient of trace molecule in medium (unit: m² s⁻¹), given
- `T` Trace medium temperature in `[K]`
- `mol` Trace molecule
- `med` Diffusion medium
"""
diffusive_coefficient(
            T::FT,
            mol::TraceGasCO₂,
            med::TraceGasAir
) where {FT<:AbstractFloat} =
(
    D_0::FT = D_CO2_IN_AIR;

    return D_0 * relative_diffusive_coefficient(T, mol, med)
)




"""
Relative diffusive coefficient of trace gas in medium

$(METHODLIST)

"""
function relative_diffusive_coefficient end




"""
    relative_diffusive_coefficient(
                T::FT,
                mol::AbstractTraceGas = TraceGasH₂O(),
                med::AbstractTraceGas = TraceGasAir()
    ) where {FT<:AbstractFloat}

Relative diffusive coefficient of trace gas in medium, given
- `T` Water vapor temperature in `[K]`
- `mol` Trace molecule. Optional, default is water vapor
- `med` Medium. Optional, default is air
"""
relative_diffusive_coefficient(
            T::FT,
            mol::AbstractTraceGas = TraceGasH₂O(),
            med::AbstractTraceGas = TraceGasAir()
) where {FT<:AbstractFloat} =
(
    return (T / K_25(FT)) ^ FT(1.8)
)








###############################################################################
#
# Latent heat
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
                med::TraceLiquidH₂O = TraceLiquidH₂O()
    ) where {FT<:AbstractFloat}

Latent heat of vaporization in `[J kg⁻¹]`, given
- `T` Medium temperature in `[K]`
- `med` Medium. Optional. Default is liquid water
"""
latent_heat_vapor(
            T::FT,
            med::TraceLiquidH₂O = TraceLiquidH₂O()
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
                med::TraceLiquidH₂O = TraceLiquidH₂O()
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
            med::TraceLiquidH₂O = TraceLiquidH₂O()
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
                med::TraceLiquidH₂O = TraceLiquidH₂O()
    ) where {FT<:AbstractFloat}

Saturation vapor pressure in `[Pa]`, given
- `T` Liquid water temperature in `[K]`
- `med` Medium. Optional. Default is liquid water
"""
saturation_vapor_pressure(
            T::FT,
            med::TraceLiquidH₂O = TraceLiquidH₂O()
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
                med::TraceLiquidH₂O = TraceLiquidH₂O()
    ) where {FT<:AbstractFloat}

Saturation vapor pressure, given
- `T` Liquid water temperature in `[K]`
- `Ψ` Liquid water pressure in `[Pa]`, positive/negative for convex/concave
    interface; if `Ψ` is given, [`pressure_correction`](@ref) is made
- `med` Medium. Optional. Default is liquid water
"""
saturation_vapor_pressure(
            T::FT,
            Ψ::FT,
            med::TraceLiquidH₂O = TraceLiquidH₂O()
) where {FT<:AbstractFloat} =
(
    return saturation_vapor_pressure(T, med) * pressure_correction(T, Ψ, med)
)




"""
First order derivative of saturation vapor pressure

$(METHODLIST)

"""
function saturation_vapor_pressure_slope end




"""
    saturation_vapor_pressure_slope(
                T::FT,
                med::TraceLiquidH₂O = TraceLiquidH₂O()
    ) where {FT<:AbstractFloat}

First order derivative of saturation vapor pressure, given
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
            med::TraceLiquidH₂O = TraceLiquidH₂O()
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
                med::TraceLiquidH₂O = TraceLiquidH₂O()
    ) where {FT<:AbstractFloat}

First order derivative of saturation vapor pressure, given
- `T` Liquid water temperature in `[K]`
- `Ψ` Liquid water pressure in `[Pa]`, positive/negative for convex/concave
    interface; if `Ψ` is given, [`pressure_correction`](@ref) is made
- `med` Medium. Optional. Default is liquid water
"""
saturation_vapor_pressure_slope(
            T::FT,
            Ψ::FT,
            med::TraceLiquidH₂O = TraceLiquidH₂O()
) where {FT<:AbstractFloat} =
(
    return saturation_vapor_pressure_slope(T, med) *
           pressure_correction(T, Ψ, med)
)








###############################################################################
#
# Surface tension functions
#
###############################################################################
const ST_k      = 0.2358
const ST_T_crit = 647.096
const ST_exp    = 1.256
const ST_corr   = 0.625
const ST_ref    = 0.07197220523

"""
Surface tension

$(METHODLIST)

"""
function surface_tension end




"""
    surface_tension(
                T::FT,
                med::TraceLiquidH₂O = TraceLiquidH₂O()
    ) where {FT<:AbstractFloat}

Surface tension `[N m⁻¹]` of trace liquid, given
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
            med::TraceLiquidH₂O = TraceLiquidH₂O()
) where {FT<:AbstractFloat} =
(
    _ST_corr    ::FT = ST_corr;
    _ST_exp     ::FT = ST_exp;
    _ST_k       ::FT = ST_k;
    _ST_T_crit  ::FT = ST_T_crit;
    _ST_T_r_diff::FT = 1 - T / _ST_T_crit;

    return _ST_k * _ST_T_r_diff^_ST_exp * (1 - _ST_corr*_ST_T_r_diff)
)




"""
Relative surface tension of trace liquid relative to 298.15 K.

$(METHODLIST)

"""
function relative_surface_tension end




"""
    relative_surface_tension(
                T::FT,
                med::TraceLiquidH₂O = TraceLiquidH₂O()
    ) where {FT<:AbstractFloat}

Relative surface tension of trace liquid relative to 298.15 K, given
- `T` Liquid water temperature in `[K]`
- `med` Medium. Optional. Default is liquid water
"""
relative_surface_tension(
            T::FT,
            med::TraceLiquidH₂O = TraceLiquidH₂O()
) where {FT<:AbstractFloat} =
(
    _ST_ref::FT = ST_ref;

    return surface_tension(T, med) / _ST_ref
)








###############################################################################
#
# Viscosity functions
#
###############################################################################
const VIS_A =  1.856e-14    # Pa s | Used for viscosity
const VIS_B =  4209.0       # K    | Used for viscosity
const VIS_C =  0.04527      # K⁻¹  | Used for viscosity
const VIS_D = -3.376e-5     # K⁻²  | Used for viscosity

"""
Viscosity.

$(METHODLIST)

"""
function viscosity end




"""
    viscosity(T::FT,
              med::TraceLiquidH₂O = TraceLiquidH₂O()
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
          med::TraceLiquidH₂O = TraceLiquidH₂O()
) where {FT<:AbstractFloat} =
(
    _A::FT = VIS_A;
    _B::FT = VIS_B;
    _C::FT = VIS_C;
    _D::FT = VIS_D;

    return _A * exp( _B/T + _C*T + _D*T^2 )
)




"""
Viscosity relative to 298.15 K

$(METHODLIST)

"""
function relative_viscosity end




"""
    relative_viscosity(
                T::FT,
                med::TraceLiquidH₂O = TraceLiquidH₂O()
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
            med::TraceLiquidH₂O = TraceLiquidH₂O()
) where {FT<:AbstractFloat} =
(
    _B::FT = VIS_B;
    _C::FT = VIS_C;
    _D::FT = VIS_D;
    _K::FT = K_25(FT);

    return exp( _B * ( 1/T - 1/_K) + _C * (T - _K) + _D * (T^2 - _K^2) )
)



end # module

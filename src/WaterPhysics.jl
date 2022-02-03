module WaterPhysics

using DocStringExtensions: METHODLIST, TYPEDEF, TYPEDFIELDS
using PkgUtility: CP_L, CP_V, GAS_R, K_BOLTZMANN, LH_V0, PRESS_TRIPLE, R_V,
      T_25, T_TRIPLE, V_H₂O
using UnPack: @unpack




# export public types
export TraceGasAir, TraceGasCO₂, TraceGasH₂O, TraceGasN₂, TraceGasO₂,
       TraceLiquidH₂O

# export public functions
export capillary_pressure, diffusive_coefficient, latent_heat_vapor,
       pressure_correction, relative_diffusive_coefficient,
       relative_surface_tension, relative_viscosity, saturation_vapor_pressure,
       saturation_vapor_pressure_slope, surface_tension, viscosity




###############################################################################
#
# Trace gas or liquid
#
###############################################################################
"""
`WaterPhysics` uses the multiple dispatch approach to calculate the temperature
    and pressure dependent physical properties of water and other molecules,
    such as CO₂. The trace molecules and mediums are catergorized to gas or
    liquid subject to a general type `AbstractTrace`. Hierarchy of
    `AbstractTrace`:

- [`AbstractTraceGas`](@ref)
- [`AbstractTraceLiquid`](@ref)

$(TYPEDEF)

"""
abstract type AbstractTrace{FT<:AbstractFloat} end




"""
The gas can be either the target trace molecule (e.g., when computing diffusive
    coefficient of CO₂ in water using [`diffusive_coefficient`](@ref)) or the
    medium (e.g., when computing diffusive coefficient or CO₂ in air using
    [`diffusive_coefficient`](@ref)). Currently, `WaterPhysics` supports the
    following subtypes of `AbstractTraceGas`:

- [`TraceGasAir`](@ref)
- [`TraceGasCO₂`](@ref)
- [`TraceGasH₂O`](@ref)
- [`TraceGasO₂`](@ref)

$(TYPEDEF)

"""
abstract type AbstractTraceGas{FT<:AbstractFloat} <: AbstractTrace{FT} end




"""
The liquid can be either the medium for gas (e.g., when computing diffusive
    coefficient of CO₂ in water using [`diffusive_coefficient`](@ref)) or the
    target molecule (e.g., when computing surface tension of water using
    [`surface_tension`](@ref)). Currently. `WaterPhysics` supports the
    following subtypes of `AbstractTraceLiquid`:

- [`TraceLiquidH₂O`](@ref)

$(TYPEDEF)

"""
abstract type AbstractTraceLiquid{FT<:AbstractFloat} <: AbstractTrace{FT} end




"""
Identity trace label for air.

$(TYPEDEF)

---
Examples
```julia
# load package before using the public functions
using WaterPhysics
# gas-phase air trace
trace = TraceGasAir{Float64}();
```
"""
struct TraceGasAir{FT<:AbstractFloat} <: AbstractTraceGas{FT} end




"""
Identity trace label for gas phase CO₂.

$(TYPEDEF)

# Fields

$(TYPEDFIELDS)

---
Examples
```julia
# gas-phase CO₂ trace
trace = TraceGasCO₂{Float64}();
```
"""
Base.@kwdef struct TraceGasCO₂{FT<:AbstractFloat} <: AbstractTraceGas{FT}
    # related to diffusive coefficient in air
    "Diffusive coefficient in air in `[m² s⁻¹]`"
    d_air  ::FT = 2.82e-5 / 1.6

    # related to diffusive coefficient in liquid water
    "Hydrodynamic radius of the solute in `[m]`"
    a_298  ::FT = 1.68e-10
    "Coefficient to make temperature correction over ydrodynamic radius"
    a_a    ::FT = 2e-3
    "Diffusive coefficient in liquid water in `[m² s⁻¹]`"
    d_water::FT = 2.147813e-9
end




"""
Identity trace label for gas phase H₂O.

$(TYPEDEF)

# Fields

$(TYPEDFIELDS)

---
Examples
```julia
# gas-phase H₂O trace
trace = TraceGasH₂O{Float64}();
```
"""
Base.@kwdef struct TraceGasH₂O{FT<:AbstractFloat} <: AbstractTraceGas{FT}
    # related to diffusive coefficient in air
    "Diffusive coefficient in air in `[m² s⁻¹]`"
    d_air::FT = 2.82e-5
end




"""
Identity trace label for gas phase N₂.

$(TYPEDEF)

# Fields

$(TYPEDFIELDS)

---
Examples
```julia
# gas-phase N₂ trace
trace = TraceGasN₂{Float64}();
```
"""
Base.@kwdef struct TraceGasN₂{FT<:AbstractFloat} <: AbstractTraceGas{FT}
    # related to diffusive coefficient in liquid water
    "Hydrodynamic radius of the solute in `[m]`"
    a_298  ::FT = 1.90e-10
    "Coefficient to make temperature correction over ydrodynamic radius"
    a_a    ::FT = 2.2e-3
    "Diffusive coefficient in liquid water in `[m² s⁻¹]`"
    d_water::FT = 1.899062e-9
end




"""
Identity trace label for gas phase O₂.

$(TYPEDEF)

# Fields

$(TYPEDFIELDS)

---
Examples
```julia
# gas-phase O₂ trace
trace = TraceGasO₂{Float64}();
```
"""
Base.@kwdef struct TraceGasO₂{FT<:AbstractFloat} <: AbstractTraceGas{FT}
    # related to diffusive coefficient in air
    "Diffusive coefficient in air in `[m² s⁻¹]`"
    d_air  ::FT = 1.76e-5

    # related to diffusive coefficient in liquid water
    "Diffusive coefficient in liquid water in `[m² s⁻¹]`"
    d_water::FT = 2.10e-9
end




"""
Identity trace label for liquid phase H₂O.

$(TYPEDEF)

# Fields

$(TYPEDFIELDS)

---
Examples
```julia
# liquid-phase water trace
trace = TraceLiquidH₂O{Float64}();
```
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
Capillary pressure of the trace liquid in a pipe is a function of surface
    tension (`γ`), pipe raduis (`r`), and contact angle (`α`):

```math
P_{c} = \\dfrac{2 ⋅ γ ⋅ \\text{cos}(α)}{r}
```

`WaterPhysics` supports the following methods to compute capillary pressure:

$(METHODLIST)

"""
function capillary_pressure end




"""
When the contact angle is not given, capillary pressure is computed using

    capillary_pressure(
                r::FT,
                T::FT,
                med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
    ) where {FT<:AbstractFloat}

Capillary pressure of trace liquid in `[Pa]`, given
- `r` Curvature radius in `[m]`
- `T` Trace liquid temperature in `[K]`
- `med` Medium. Optional. Default is liquid water

---
Examples
```julia
# capillary pressure of liquid water
p_1 = capillary_pressure(3e-6, 298.15);
p_2 = capillary_pressure(3e-6, 298.15, TraceLiquidH₂O{Float64}());
```
"""
capillary_pressure(
            r::FT,
            T::FT,
            med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
    2 * surface_tension(T, med) / r




"""
When the contact angle is given, capillary pressure is computed using

    capillary_pressure(
                r::FT,
                T::FT,
                α::FT,
                med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
    ) where {FT<:AbstractFloat}

Capillary pressure of trace liquid in `[Pa]`, given
- `r` Pipe radius in `[m]`
- `T` Trace liquid temperature in `[K]`
- `α` Contact angle in `[°]`
- `med` Medium. Optional. Default is liquid water

---
Examples
```julia
# capillary pressure of liquid water with contact angle
p_1 = capillary_pressure(3e-6, 298.15, 30.0);
p_2 = capillary_pressure(3e-6, 298.15, 30.0, TraceLiquidH₂O{Float64}());
```
"""
capillary_pressure(
            r::FT,
            T::FT,
            α::FT,
            med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
    cosd(α) * capillary_pressure(r, T, med)








###############################################################################
#
# Diffusive coefficient
#
###############################################################################
"""
Diffusion of trace molecules in medium is temperature dependent, to calculate
    this temperature dependency, we provide a function to quick estimate this
    value for different trace molecules using `diffusive_coefficient`. The
    supported methods are:

$(METHODLIST)

"""
function diffusive_coefficient end




"""
For the diffusive coefficient of gas in air, the coefficient is simply treated
    as a function of temperature and the reference coefficient at 298.15 K:

```@math
D = D_\\text{ref} ⋅ \\left( \\dfrac{T}{298.15} \\right)^{1.8}
```

This temperature correction is also available as a stand alone function:
    [`relative_diffusive_coefficient`](@ref).

    diffusive_coefficient(
                T::FT,
                mol::Union{TraceGasCO₂{FT}, TraceGasH₂O{FT}, TraceGasO₂{FT}},
                med::TraceGasAir{FT}
    ) where {FT<:AbstractFloat}

Diffusive coefficient of trace molecule in `[m² s⁻¹]`, given
- `T` Trace medium temperature in `[K]`
- `mol` Trace molecule
- `med` Diffusion medium (air)

---
Examples
```julia
# diffusive coefficient of gas-phase CO₂, H₂O, and O₂ in air
air   = TraceGasAir{Float64}();
d_CO₂ = diffusive_coefficient(298.15, TraceGasCO₂{Float64}(), air);
d_H₂O = diffusive_coefficient(298.15, TraceGasH₂O{Float64}(), air);
d_O₂  = diffusive_coefficient(298.15, TraceGasO₂{Float64}(), air);
```
"""
diffusive_coefficient(
            T::FT,
            mol::Union{TraceGasCO₂{FT}, TraceGasH₂O{FT}, TraceGasO₂{FT}},
            med::TraceGasAir{FT}
) where {FT<:AbstractFloat} =
    mol.d_air * relative_diffusive_coefficient(T, mol, med)




"""
The diffusive coefficient of disolved gas in liquid water is also a function of
    temperature. The calculation is based on the empirical formulation and
    variables from Cadogen et al. (2015) Diffusion coefficients of CO₂ and N₂
    in water at temperatures between 298.15 K and 423.15 K at pressures up to
    45 MPa.

    diffusive_coefficient(
                T::FT,
                mol::Union{TraceGasCO₂{FT}, TraceGasN₂{FT}},
                med::TraceLiquidH₂O{FT}
    ) where {FT<:AbstractFloat}

Diffusive coefficient of trace molecule in `[m² s⁻¹]`, given
- `T` Trace medium temperature in `[K]`
- `mol` Trace molecule
- `med` Diffusion medium (liquid water)

---
Examples
```julia
# diffusive coefficient of CO₂ and N₂ in liquid water
water = TraceLiquidH₂O{Float64}();
d_CO₂ = diffusive_coefficient(298.15, TraceGasCO₂{FT}(), water);
d_N₂  = diffusive_coefficient(298.15, TraceGasN₂{FT}(), water);
```
"""
diffusive_coefficient(
            T::FT,
            mol::Union{TraceGasCO₂{FT}, TraceGasN₂{FT}},
            med::TraceLiquidH₂O{FT}
) where {FT<:AbstractFloat} =
(
    @unpack a_298, a_a = mol;
    _a = a_298 * (1 + a_a *  (T-298));

    return K_BOLTZMANN(FT) * T / (4 * FT(pi) * viscosity(T, med) * _a)
)




"""
The temperature dependency of diffusive coefficient impacts leaf gas exchange
    via change the maximal stomatal conductance to H₂O vapor and CO₂, given
    that Since the stomatal conductance should not exceed its structural limit.
    To account for this effect, we provide a function to calculate the
    diffusive coefficient relative to 25 Celcius.

$(METHODLIST)

"""
function relative_diffusive_coefficient end




"""
As mentioned in [`diffusive_coefficient`](@ref), relative diffusive coefficient
    of target gas in gas medium is calculated as an exponential function of
    temperature. The shortcut method is

    relative_diffusive_coefficient(
                T::FT,
                mol::AbstractTraceGas{FT} = TraceGasH₂O{FT}(),
                med::AbstractTraceGas{FT} = TraceGasAir{FT}()
    ) where {FT<:AbstractFloat}

Relative diffusive coefficient of trace gas in medium, given
- `T` Water vapor temperature in `[K]`
- `mol` Trace molecule. Optional, default is water vapor
- `med` Medium. Optional, default is air

---
Examples
```julia
# relative diffusive coefficient of gas in air
rat_1 = relative_diffusive_coefficient(298.15);
rat_2 = relative_diffusive_coefficient(298.15, TraceGasH₂O{Float64}());
rat_3 = relative_diffusive_coefficient(298.15, TraceGasH₂O{Float64}(),
                                       TraceGasAir{Float64});
```
"""
relative_diffusive_coefficient(
            T::FT,
            mol::AbstractTraceGas{FT} = TraceGasH₂O{FT}(),
            med::AbstractTraceGas{FT} = TraceGasAir{FT}()
) where {FT<:AbstractFloat} =
    (T / T_25(FT)) ^ FT(1.8)




"""
As to disolved gas diffsuin in liquid medium, the relative coefficient is also
    computed based on the empirical formulation from Cadogan et al. (2014).

    relative_diffusive_coefficient(
                T::FT,
                mol::AbstractTraceGas{FT}
                med::AbstractTraceGas{FT}
    ) where {FT<:AbstractFloat}

Relative diffusive coefficient of trace gas in medium, given
- `T` Water vapor temperature in `[K]`
- `mol` Trace molecule. Optional, default is water vapor
- `med` Medium. Optional, default is air

---
Examples
```julia
# relative diffusive coefficient of CO₂ and N₂ in liquid water
water = TraceLiquidH₂O{Float64}();
rat_1 = relative_diffusive_coefficient(298.15, TraceGasCO₂{Float64}(), water);
rat_2 = relative_diffusive_coefficient(298.15, TraceGasN₂{Float64}(), water);
```
"""
relative_diffusive_coefficient(
            T::FT,
            mol::Union{TraceGasCO₂{FT}, TraceGasN₂{FT}},
            med::TraceLiquidH₂O{FT}
) where {FT<:AbstractFloat} =
    diffusive_coefficient(T, mol, med) / mol.d_water








###############################################################################
#
# Latent heat of evaporation
# Adapted from the ClimateMachine/Common/ThermoDynamics/relations.jl
#
###############################################################################
"""
Water evaporation from liquid phase is a key process to regulate leaf
    temperature, and to best represent this process. We computed the latent
    heat coefficient from water temperature:

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

---
Examples
```julia
# latent heat of vapor of liquid water
λ_1 = latent_heat_vapor(298.15);
λ_2 = latent_heat_vapor(298.15, TraceLiquidH₂O{Float64}());
```
"""
latent_heat_vapor(
            T::FT,
            med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
    LH_V0(FT) + (CP_V(FT) - CP_L(FT)) * (T - T_TRIPLE(FT))








###############################################################################
#
# Saturated vapor pressure
# Adapted from the ClimateMachine/Common/ThermoDynamics/relations.jl
#
###############################################################################
"""
Yet, the saturation vapor pressure is not only a function of temperature, but
    also a function of the air-water interface curvature, known as the Kelvin
    equation. The package provide [`pressure_correction`](@ref) to make the
    correction.

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

---
Examples
```julia
# correction to saturation vapor pressure of liquid water with non-flat surface
cor_1 = pressure_correction(298.15, -1e6);
cor_2 = pressure_correction(298.15, -1e6, TraceLiquidH₂O{Float64}());
```
"""
pressure_correction(
            T::FT,
            Ψ::FT,
            med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
    exp((Ψ * V_H₂O(FT)) / (GAS_R(FT) * T))




"""
When temperature increases, liquid water vapor pressure increases
    exponentially. And this correlation is accounted for using the functions
    below:

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

---
Examples
```julia
# saturation vapor pressure of liquid water with flat surface
p_1 = saturation_vapor_pressure(298.15);
p_2 = saturation_vapor_pressure(298.15, TraceLiquidH₂O{Float64}());
```
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

---
Examples
```julia
# saturation vapor pressure of liquid water with non-flat surface
p_1 = saturation_vapor_pressure(298.15, -1e6);
p_2 = saturation_vapor_pressure(298.15, -1e6, TraceLiquidH₂O{Float64}());
```
"""
saturation_vapor_pressure(
            T::FT,
            Ψ::FT,
            med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
    saturation_vapor_pressure(T, med) * pressure_correction(T, Ψ, med)




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

---
Examples
```julia
# slope of saturation vapor pressure of liquid water with flat surface
s_1 = saturation_vapor_pressure_slope(298.15);
s_2 = saturation_vapor_pressure_slope(298.15, TraceLiquidH₂O{Float64}());
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

---
Examples
```julia
# slope of saturation vapor pressure of liquid water with non-flat surface
s_1 = saturation_vapor_pressure_slope(298.15, -1e6);
s_2 = saturation_vapor_pressure_slope(298.15, -1e6, TraceLiquidH₂O{Float64}());
```
"""
saturation_vapor_pressure_slope(
            T::FT,
            Ψ::FT,
            med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
    saturation_vapor_pressure_slope(T, med) * pressure_correction(T, Ψ, med)








###############################################################################
#
# Surface tension
#
###############################################################################
"""
When water temperature increases, the surface tension at the air-water
    interface decreases. Surface tension changes impacts the plant water
    transport via two aspects. First, if surface tension is lower, for a
    constant soil water content, the soil matrix potential gets less negative
    because the capillary pressure at the air-water interface decreases. And
    this is beneficial to plants. Second, the air-water interface at the pit
    membrane also has lower capillary pressure when temperature increases,
    meaning that the xylem conduits are less resistant to drought-induced
    air-seeded cavitation. And this is harmful for plants. Though the surface
    tension does not differ much with temperature change within the plant
    physiological active range, we account for this effect in our Land model.

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

---
Examples
```julia
# surface tension of liquid water
γ_1 = surface_tension(298.15);
γ_2 = surface_tension(298.15, TraceLiquidH₂O{Float64}());
```
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

---
Examples
```julia
# relative surface tension of liquid water
γ_1 = relative_surface_tension(298.15);
γ_2 = relative_surface_tension(298.15, TraceLiquidH₂O{Float64}());
```
"""
relative_surface_tension(
            T::FT,
            med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
    surface_tension(T, med) / med.γ_ref








###############################################################################
#
# Viscosity
#
###############################################################################
"""
When temperature increases, liquid water viscosuty decreases, meaning that the
    resistance for water decreases and the pressure drop per flow rate
    decreases. This effect is pretty significant as 1 degree increase of
    temperature results in about 2.4% drop in viscosity, and this is very
    beneficial to plant water transport. Unfortunately, to our knowledge, very
    few models account for this effect when modeling plant hydraulics because
    of the difficulty in modeling the energy budget along the flow path. We
    plan to have this effect accounted for in our CliMA Land model, by
    computing the water tempreature along the flow path and thus the viscosity
    change. The functions to be used are:

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

---
Examples
```julia
# viscosity of liquid water
γ_1 = viscosity(298.15);
γ_2 = viscosity(298.15, TraceLiquidH₂O{Float64}());
```
"""
viscosity(T::FT,
          med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
    med.υ_A * exp( med.υ_B/T + med.υ_C*T + med.υ_D*T^2 )




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

---
Examples
```julia
# relative viscosity of liquid water
γ_1 = relative_viscosity(298.15);
γ_2 = relative_viscosity(298.15, TraceLiquidH₂O{Float64}());
```
"""
relative_viscosity(
            T::FT,
            med::TraceLiquidH₂O{FT} = TraceLiquidH₂O{FT}()
) where {FT<:AbstractFloat} =
(
    @unpack υ_B, υ_C, υ_D = med;
    _K = T_25(FT);

    return exp( υ_B * ( 1/T - 1/_K) + υ_C * (T - _K) + υ_D * (T^2 - _K^2) )
)




end # module

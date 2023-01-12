#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-14: Move the structure from Photosynthesis.jl, only P_A and P_O2 for now
#     2022-Jan-14: rename P_A to P_AIR, P_O2 to P_O₂
#     2022-Jan-14: add p_CO₂ to the structure
#     2022-Mar-09: add t, p_H₂O, p_H₂O_sat, rh, and wind fields
#     2022-Jul-13: remove fields p_H₂O_sat and rh to avoid update issues
#     2022-Jul-20: remove fields P_O₂ to avoid update issues
#     2022-Jul-20: add fields: Z, ΔZ, e, n_CO₂, n_H₂O, ∂e∂t, ∂CO₂∂t, and ∂H₂O∂t
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores air layer information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct AirLayer{FT<:AbstractFloat}
    # Location and geometry of the air layer
    "Mean height of the layer `[m]`"
    Z::FT = 0.5
    "Layer thickness `[m]`"
    ΔZ::FT = 1

    # Parameters that are not supposed to change with time
    "Atmospheric pressure `[Pa]`"
    P_AIR::FT = P_ATM()

    # Prognostic variables (not used for ∂y∂t)
    "CO₂ partial pressure `[Pa]`"
    p_CO₂::FT = 40
    "H₂O partial pressure `[Pa]`"
    p_H₂O::FT = 1500
    "Temperature `[K]`"
    t::FT = T₂₅()
    "Wind speed `[m s⁻¹]`"
    wind::FT = 1

    # Prognostic variables (used for ∂y∂t)
    "Total energy within the air layer `[J m⁻²]`"
    e::FT = CP_D_MOL() * (P_AIR - p_H₂O) * ΔZ / GAS_R() + CP_V_MOL() * p_H₂O * ΔZ / GAS_R()
    "Mole of CO₂ per surface area `[mol m⁻²]`"
    n_CO₂::FT = p_CO₂ * ΔZ / (GAS_R() * t)
    "Mole of H₂O per surface area `[mol m⁻²]`"
    n_H₂O::FT = p_H₂O * ΔZ / (GAS_R() * t)
    "Marginal increase in total energy `[J m⁻² s⁻¹]`"
    ∂e∂t::FT = 0
    "Marginal increase in total moles of CO₂ `[mol m⁻² s⁻¹]`"
    ∂CO₂∂t::FT = 0
    "Marginal increase in total moles of H₂O `[mol m⁻² s⁻¹]`"
    ∂H₂O∂t::FT = 0
end

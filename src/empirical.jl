#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-01: migrate function from older version
#     2022-Jul-01: rename the function from stomatal_conductance to empirical_equation
#
#######################################################################################################################################################################################################
"""
This function returns the stomatal conductance computed from empirical stomatal models. This is not the solution! Supported methods are

$(METHODLIST)

"""
function empirical_equation end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add method for BallBerrySM
#
#######################################################################################################################################################################################################
"""

    empirical_equation(sm::BallBerrySM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat}

Return the stomatal conductance computed from empirical model formulation, given
- `sm` `BallBerrySM` type model
- `leaf` `Leaf` type struct
- `air` `AirLayer` type environmental conditions
- `β` Tuning factor for G1 (must be 1 if tuning factor is not based on G1)
"""
empirical_equation(sm::BallBerrySM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack G0, G1 = sm;
    @unpack PSM = leaf;
    @unpack P_AIR = air;

    return G0 + β * G1 * air.rh * PSM.a_net * FT(1e-6) / leaf.p_CO₂_s * P_AIR
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add method for GentineSM
#
#######################################################################################################################################################################################################
"""

    empirical_equation(sm::GentineSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat}

Return the stomatal conductance computed from empirical model formulation, given
- `sm` `GentineSM` type model
- `leaf` `Leaf` type struct
- `air` `AirLayer` type environmental conditions
- `β` Tuning factor for G1 (must be 1 if tuning factor is not based on G1)
"""
empirical_equation(sm::GentineSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack G0, G1 = sm;
    @unpack P_AIR = air;

    return G0 + β * G1 * leaf.PSM.a_net * FT(1e-6) / leaf.p_CO₂_i * P_AIR
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add method for LeuningSM
#
#######################################################################################################################################################################################################
"""

    empirical_equation(sm::LeuningSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat}

Return the stomatal conductance computed from empirical model formulation, given
- `sm` `LeuningSM` type model
- `leaf` `Leaf` type struct
- `air` `AirLayer` type environmental conditions
- `β` Tuning factor for G1 (must be 1 if tuning factor is not based on G1)
"""
empirical_equation(sm::LeuningSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack D0, G0, G1 = sm;
    @unpack PSM = leaf;
    @unpack P_AIR = air;

    _γ_s = (typeof(PSM) <: C4VJPModel) ? 0 : PSM.γ_star;

    return G0 + β * G1 / (1 + (leaf.p_H₂O_sat - air.p_H₂O) / D0) * PSM.a_net * FT(1e-6) / (leaf.p_CO₂_s - _γ_s) * P_AIR
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add method for MedlynSM
# Bug fix:
#     2022-Jul-07: add the factor 1.6 for Medlyn model
#
#######################################################################################################################################################################################################
"""

    empirical_equation(sm::MedlynSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat}

Return the stomatal conductance computed from empirical model formulation, given
- `sm` `MedlynSM` type model
- `leaf` `Leaf` type struct
- `air` `AirLayer` type environmental conditions
- `β` Tuning factor for G1 (must be 1 if tuning factor is not based on G1)
"""
empirical_equation(sm::MedlynSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack G0, G1 = sm;
    @unpack PSM = leaf;
    @unpack P_AIR = air;

    return G0 + FT(1.6) * (1 + β * G1 / sqrt(leaf.p_H₂O_sat - air.p_H₂O)) * leaf.PSM.a_net * FT(1e-6) / air.p_CO₂ * P_AIR
);

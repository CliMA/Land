#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-01: migrate function from older version
#     2022-Jul-01: rename the function from stomatal_conductance to empirical_equation
#
#######################################################################################################################################################################################################
"""
This function returns the stomatal conductance computed from empirical stomatal models. This is not the solution! Supported methods are for
- Leaf
- Leaves1D (ind=1 for sunlit, ind=2 for shaded leaves)
- Leaves2D (ind=NA for shaded, ind>1 for sunlit leaves)

"""
function empirical_equation end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add method for BallBerrySM
#     2022-Jul-07: use a_net stored in Leaf
#     2022-Jul-07: add method for GentineSM
#     2022-Jul-07: add method for LeuningSM
#     2022-Jul-07: add method for MedlynSM
#     2022-Oct-20: add a max controller to make sure vpd is at least 1 Pa
# Bug fix:
#     2022-Jul-07: add the factor 1.6 for Medlyn model
#
#######################################################################################################################################################################################################
"""

    empirical_equation(sm::BallBerrySM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat}
    empirical_equation(sm::GentineSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat}
    empirical_equation(sm::LeuningSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat}
    empirical_equation(sm::MedlynSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat}

Return the stomatal conductance computed from empirical model formulation, given
- `sm` `BallBerrySM`, `GentineSM`, `LeuningSM`, or `MedlynSM` type model
- `leaf` `Leaf` type struct
- `air` `AirLayer` type environmental conditions
- `β` Tuning factor for G1 (must be 1 if tuning factor is not based on G1)

"""
empirical_equation(sm::BallBerrySM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack G0, G1 = sm;
    @unpack P_AIR = air;

    return G0 + β * G1 * air.p_H₂O / saturation_vapor_pressure(air.t) * leaf.a_net * FT(1e-6) / leaf._p_CO₂_s * P_AIR
);

empirical_equation(sm::GentineSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack G0, G1 = sm;
    @unpack P_AIR = air;

    return G0 + β * G1 * leaf.a_net * FT(1e-6) / leaf._p_CO₂_i * P_AIR
);

empirical_equation(sm::LeuningSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack D0, G0, G1 = sm;
    @unpack PSM = leaf;
    @unpack P_AIR = air;

    _γ_s = (typeof(PSM) <: C4VJPModel) ? 0 : PSM._γ_star;
    _vpd = max(1, saturation_vapor_pressure(leaf.t) - air.p_H₂O);

    return G0 + β * G1 / (1 + _vpd / D0) * leaf.a_net * FT(1e-6) / (leaf._p_CO₂_s - _γ_s) * P_AIR
);

empirical_equation(sm::MedlynSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack G0, G1 = sm;
    @unpack P_AIR = air;

    _vpd = max(1, saturation_vapor_pressure(leaf.t) - air.p_H₂O);

    return G0 + FT(1.6) * (1 + β * G1 / sqrt(_vpd)) * leaf.a_net * FT(1e-6) / air.p_CO₂ * P_AIR
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add method for BallBerrySM using Leaves1D
#     2022-Jul-07: add method for GentineSM using Leaves1D
#     2022-Jul-07: add method for LeuningSM using Leaves1D
#     2022-Jul-07: add method for MedlynSM using Leaves1D
#     2022-Oct-20: add a max controller to make sure vpd is at least 1 Pa
#
#######################################################################################################################################################################################################
"""

    empirical_equation(sm::BallBerrySM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat}
    empirical_equation(sm::GentineSM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat}
    empirical_equation(sm::LeuningSM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat}
    empirical_equation(sm::MedlynSM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat}

Return the stomatal conductance computed from empirical model formulation, given
- `sm` `BallBerrySM`, `GentineSM`, `LeuningSM`, or `MedlynSM` type model
- `leaves` `Leaves1D` type struct
- `air` `AirLayer` type environmental conditions
- `ind` Leaf index (1 for sunlit and 2 for shaded)
- `β` Tuning factor for G1 (must be 1 if tuning factor is not based on G1)

"""
empirical_equation(sm::BallBerrySM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack G0, G1 = sm;
    @unpack P_AIR = air;

    return G0 + β * G1 * air.p_H₂O / saturation_vapor_pressure(air.t) * leaves.a_net[ind] * FT(1e-6) / leaves._p_CO₂_s[ind] * P_AIR
);

empirical_equation(sm::GentineSM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack G0, G1 = sm;
    @unpack P_AIR = air;

    return G0 + β * G1 * leaves.a_net[ind] * FT(1e-6) / leaves._p_CO₂_i[ind] * P_AIR
);

empirical_equation(sm::LeuningSM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack D0, G0, G1 = sm;
    @unpack P_AIR = air;

    _γ_s = (typeof(leaves.PSM) <: C4VJPModel) ? 0 : leaves.PSM._γ_star;
    _vpd = max(1, saturation_vapor_pressure(leaves.t[ind]) - air.p_H₂O);

    return G0 + β * G1 / (1 + _vpd / D0) * leaves.a_net[ind] * FT(1e-6) / (leaves._p_CO₂_s[ind] - _γ_s) * P_AIR
);

empirical_equation(sm::MedlynSM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack G0, G1 = sm;
    @unpack P_AIR = air;

    _vpd = max(1, saturation_vapor_pressure(leaves.t[ind]) - air.p_H₂O);

    return G0 + FT(1.6) * (1 + β * G1 / sqrt(_vpd)) * leaves.a_net[ind] * FT(1e-6) / air.p_CO₂ * P_AIR
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add method for BallBerrySM using Leaves2D for shaded leaves
#     2022-Jul-07: add method for GentineSM using Leaves2D for shaded leaves
#     2022-Jul-07: add method for LeuningSM using Leaves2D for shaded leaves
#     2022-Jul-07: add method for MedlynSM using Leaves2D for shaded leaves
#     2022-Oct-20: add a max controller to make sure vpd is at least 1 Pa
#
#######################################################################################################################################################################################################
"""

    empirical_equation(sm::BallBerrySM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat}
    empirical_equation(sm::GentineSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat}
    empirical_equation(sm::LeuningSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat}
    empirical_equation(sm::MedlynSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat}

Return the stomatal conductance computed from empirical model formulation for the shaded leaves of `Leaves2D`, given
- `sm` `BallBerrySM`, `GentineSM`, `LeuningSM`, or `MedlynSM` type model
- `leaves` `Leaves2D` type struct
- `air` `AirLayer` type environmental conditions
- `β` Tuning factor for G1 (must be 1 if tuning factor is not based on G1)

"""
empirical_equation(sm::BallBerrySM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack G0, G1 = sm;
    @unpack P_AIR = air;

    return G0 + β * G1 * air.p_H₂O / saturation_vapor_pressure(air.t) * leaves.a_net_shaded * FT(1e-6) / leaves._p_CO₂_s_shaded * P_AIR
);

empirical_equation(sm::GentineSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack G0, G1 = sm;
    @unpack P_AIR = air;

    return G0 + β * G1 * leaves.a_net_shaded * FT(1e-6) / leaves._p_CO₂_i_shaded * P_AIR
);

empirical_equation(sm::LeuningSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack D0, G0, G1 = sm;
    @unpack P_AIR = air;

    _γ_s = (typeof(leaves.PSM) <: C4VJPModel) ? 0 : leaves.PSM._γ_star;
    _vpd = max(1, saturation_vapor_pressure(leaves.t) - air.p_H₂O);

    return G0 + β * G1 / (1 + _vpd / D0) * leaves.a_net_shaded * FT(1e-6) / (leaves._p_CO₂_s_shaded - _γ_s) * P_AIR
);

empirical_equation(sm::MedlynSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack G0, G1 = sm;
    @unpack P_AIR = air;

    _vpd = max(1, saturation_vapor_pressure(leaves.t) - air.p_H₂O);

    return G0 + FT(1.6) * (1 + β * G1 / sqrt(_vpd)) * leaves.a_net_shaded * FT(1e-6) / air.p_CO₂ * P_AIR
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add method for BallBerrySM using Leaves2D for sunlit leaves
#     2022-Jul-07: add method for GentineSM using Leaves2D for sunlit leaves
#     2022-Jul-07: add method for LeuningSM using Leaves2D for sunlit leaves
#     2022-Jul-07: add method for MedlynSM using Leaves2D for sunlit leaves
#     2022-Oct-20: add a max controller to make sure vpd is at least 1 Pa
#
#######################################################################################################################################################################################################
"""

    empirical_equation(sm::BallBerrySM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat}
    empirical_equation(sm::GentineSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat}
    empirical_equation(sm::LeuningSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat}
    empirical_equation(sm::MedlynSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat}

Return the stomatal conductance computed from empirical model formulation for the sunlit leaves of `Leaves2D`, given
- `sm` `BallBerrySM`, `GentineSM`, `LeuningSM`, or `MedlynSM` type model
- `leaves` `Leaves2D` type struct
- `air` `AirLayer` type environmental conditions
- `ind` Sunlit leaf index within the leaf angular distribution
- `β` Tuning factor for G1 (must be 1 if tuning factor is not based on G1)

"""
empirical_equation(sm::BallBerrySM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack G0, G1 = sm;
    @unpack P_AIR = air;

    return G0 + β * G1 * air.p_H₂O / saturation_vapor_pressure(air.t) * leaves.a_net_sunlit[ind] * FT(1e-6) / leaves._p_CO₂_s_sunlit[ind] * P_AIR
);

empirical_equation(sm::GentineSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack G0, G1 = sm;
    @unpack P_AIR = air;

    return G0 + β * G1 * leaves.a_net_sunlit[ind] * FT(1e-6) / leaves._p_CO₂_i_sunlit[ind] * P_AIR
);

empirical_equation(sm::LeuningSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack D0, G0, G1 = sm;
    @unpack P_AIR = air;

    _γ_s = (typeof(leaves.PSM) <: C4VJPModel) ? 0 : leaves.PSM._γ_star;
    _vpd = max(1, saturation_vapor_pressure(leaves.t) - air.p_H₂O);

    return G0 + β * G1 / (1 + _vpd / D0) * leaves.a_net_sunlit[ind] * FT(1e-6) / (leaves._p_CO₂_s_sunlit[ind] - _γ_s) * P_AIR
);

empirical_equation(sm::MedlynSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack G0, G1 = sm;
    @unpack P_AIR = air;

    _vpd = max(1, saturation_vapor_pressure(leaves.t) - air.p_H₂O);

    return G0 + FT(1.6) * (1 + β * G1 / sqrt(_vpd)) * leaves.a_net_sunlit[ind] * FT(1e-6) / air.p_CO₂ * P_AIR
);

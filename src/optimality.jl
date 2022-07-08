#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-07: migrate function from older version
# To do
#     TODO: figure out why the ratios are 1.35 and 1.6, and make them more accurate
#
#######################################################################################################################################################################################################
"""
This function returns the partial derivative of photosynthetic rate per transpiration rate. Supported methods are

$(METHODLIST)

"""
function ∂A∂E end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add the method for Leaf and Leaves2D (shaded leaves)
#
#######################################################################################################################################################################################################
"""

    ∂A∂E(leaf::Leaf{FT}, air::AirLayer{FT}) where {FT<:AbstractFloat}
    ∂A∂E(leaves::Leaves2D{FT}, air::AirLayer{FT}) where {FT<:AbstractFloat}

Return the partial derivative of A per E, given
- `leaf` `Leaf`, and `Leaves2D` (shaded fraction) type leaf
- `air` `AirLayer` type environmental conditions
"""
∂A∂E(leaf::Leaf{FT}, air::AirLayer{FT}) where {FT<:AbstractFloat} = (
    @unpack P_AIR = air;

    # compute the A and E at the current setting
    _gs1 = leaf.g_H₂O_s;
    _gh1 = 1 / (1 / _gs1 + 1 / (FT(1.35) * leaf.g_CO₂_b));
    _e1  = _gh1 * (leaf.p_H₂O_sat - air.p_H₂O) / P_AIR;
    _a1  = leaf.a_net;

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    _gs2 = _gs1 + FT(0.0001);
    _gh2 = 1 / (1 / _gs2 + 1 / (FT(1.35) * leaf.g_CO₂_b));
    _gc2 = 1 / (FT(1.6) / _gs2 + 1 / leaf.g_CO₂_b);
    leaf_photosynthesis!(leaf, air, _gc2, leaf.ppar);
    _e2 = _gh2 * (leaf.p_H₂O_sat - air.p_H₂O) / P_AIR;
    _a2 = leaf.PSM.a_net;

    return (_a2 - _a1) / (_e2 - _e1)
);


∂A∂E(leaves::Leaves2D{FT}, air::AirLayer{FT}) where {FT<:AbstractFloat} = (
    @unpack P_AIR = air;

    # compute the A and E at the current setting
    _gs1 = leaves.g_H₂O_s_shaded;
    _gh1 = 1 / (1 / _gs1 + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e1  = _gh1 * (leaves.p_H₂O_sat - air.p_H₂O) / P_AIR;
    _a1  = leaves.a_net_shaded;

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    _gs2 = _gs1 + FT(0.0001);
    _gh2 = 1 / (1 / _gs2 + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _gc2 = 1 / (FT(1.6) / _gs2 + 1 / leaves.g_CO₂_b);
    leaf_photosynthesis!(leaves, air, _gc2, leaves.ppar_shaded);
    _e2 = _gh2 * (leaves.p_H₂O_sat - air.p_H₂O) / P_AIR;
    _a2 = leaves.PSM.a_net;

    return (_a2 - _a1) / (_e2 - _e1)
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add the method for Leaves1D and Leaves2D (sunlit leaves)
#
#######################################################################################################################################################################################################
"""

    ∂A∂E(leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int) where {FT<:AbstractFloat}
    ∂A∂E(leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int) where {FT<:AbstractFloat}

Return the partial derivative of A per E, given
- `leaf` `Leaves1D`, and `Leaves2D` (sunlit fraction) type leaf
- `air` `AirLayer` type environmental conditions
- `ind` Index of the leaves (1 for sunlit and 2 for shaded for Leaves1D, all sunlit for Leaves2D)
"""
∂A∂E(leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int) where {FT<:AbstractFloat} = (
    @unpack P_AIR = air;

    # compute the A and E at the current setting
    _gs1 = leaves.g_H₂O_s[ind];
    _gh1 = 1 / (1 / _gs1 + 1 / (FT(1.35) * leaves.g_CO₂_b[ind]));
    _e1  = _gh1 * (leaves.p_H₂O_sat[ind] - air.p_H₂O) / P_AIR;
    _a1  = leaves.a_net[ind];

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    _gs2 = _gs1 + FT(0.0001);
    _gh2 = 1 / (1 / _gs2 + 1 / (FT(1.35) * leaves.g_CO₂_b[ind]));
    _gc2 = 1 / (FT(1.6) / _gs2 + 1 / leaves.g_CO₂_b[ind]);
    leaf_photosynthesis!(leaves, air, _gc2, leaves.ppar[ind]);
    _e2 = _gh2 * (leaves.p_H₂O_sat[ind] - air.p_H₂O) / P_AIR;
    _a2 = leaves.PSM.a_net;

    return (_a2 - _a1) / (_e2 - _e1)
);


∂A∂E(leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int) where {FT<:AbstractFloat} = (
    @unpack P_AIR = air;

    # compute the A and E at the current setting
    _gs1 = leaves.g_H₂O_s_sunlit[ind];
    _gh1 = 1 / (1 / _gs1 + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e1  = _gh1 * (leaves.p_H₂O_sat - air.p_H₂O) / P_AIR;
    _a1  = leaves.a_net_sunlit[ind];

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    _gs2 = _gs1 + FT(0.0001);
    _gh2 = 1 / (1 / _gs2 + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _gc2 = 1 / (FT(1.6) / _gs2 + 1 / leaves.g_CO₂_b);
    leaf_photosynthesis!(leaves, air, _gc2, leaves.ppar_sunlit[ind]);
    _e2 = _gh2 * (leaves.p_H₂O_sat - air.p_H₂O) / P_AIR;
    _a2 = leaves.PSM.a_net;

    return (_a2 - _a1) / (_e2 - _e1)
);



function ∂Θ∂E end

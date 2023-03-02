#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-07: migrate function from older version
#     2022-Jul-07: add the method for Leaf
#     2022-Jul-07: add the method for Leaves1D
#     2022-Jul-07: add the method for Leaves2D (shaded leaves)
#     2022-Jul-07: add the method for Leaves2D (sunlit leaves)
#     2022-Jul-11: deflate documentations
# To do
#     TODO: figure out why the ratios are 1.35 and 1.6, and make them more accurate
#
#######################################################################################################################################################################################################
"""

    ∂A∂E(leaf::Leaf{FT}, air::AirLayer{FT}) where {FT<:AbstractFloat}
    ∂A∂E(leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int) where {FT<:AbstractFloat}
    ∂A∂E(leaves::Leaves2D{FT}, air::AirLayer{FT}) where {FT<:AbstractFloat}
    ∂A∂E(leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int) where {FT<:AbstractFloat}

Return the partial derivative of A per E, given
- `leaf` `Leaf` type leaf
- `air` `AirLayer` type environmental conditions
- `leaves` `Leaves1D`, and `Leaves2D` type leaf
- `ind` Index of the leaves (1 for sunlit and 2 for shaded for Leaves1D, all sunlit for Leaves2D)

"""
function ∂A∂E end

∂A∂E(leaf::Leaf{FT}, air::AirLayer{FT}) where {FT<:AbstractFloat} = (
    (; P_AIR) = air;

    # compute the A and E at the current setting
    _gs1 = leaf.g_H₂O_s;
    _gh1 = 1 / (1 / _gs1 + 1 / (FT(1.35) * leaf.g_CO₂_b));
    _e1  = _gh1 * (saturation_vapor_pressure(leaf.t) - air.p_H₂O) / P_AIR;
    _a1  = leaf.a_net;

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    _gs2 = _gs1 + FT(0.0001);
    _gh2 = 1 / (1 / _gs2 + 1 / (FT(1.35) * leaf.g_CO₂_b));
    _gc2 = 1 / (FT(1.6) / _gs2 + 1 / leaf.g_CO₂_b);
    leaf_photosynthesis!(leaf, air, _gc2, leaf.ppar, leaf.t);
    _e2 = _gh2 * (saturation_vapor_pressure(leaf.t) - air.p_H₂O) / P_AIR;
    _a2 = leaf.PSM.a_net;

    return (_a2 - _a1) / (_e2 - _e1)
);

∂A∂E(leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int) where {FT<:AbstractFloat} = (
    (; P_AIR) = air;

    # compute the A and E at the current setting
    _gs1 = leaves.g_H₂O_s[ind];
    _gh1 = 1 / (1 / _gs1 + 1 / (FT(1.35) * leaves.g_CO₂_b[ind]));
    _e1  = _gh1 * (saturation_vapor_pressure(leaves.t[ind]) - air.p_H₂O) / P_AIR;
    _a1  = leaves.a_net[ind];

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    _gs2 = _gs1 + FT(0.0001);
    _gh2 = 1 / (1 / _gs2 + 1 / (FT(1.35) * leaves.g_CO₂_b[ind]));
    _gc2 = 1 / (FT(1.6) / _gs2 + 1 / leaves.g_CO₂_b[ind]);
    leaf_photosynthesis!(leaves, air, _gc2, leaves.ppar[ind], leaves.t[ind]);
    _e2 = _gh2 * (saturation_vapor_pressure(leaves.t[ind]) - air.p_H₂O) / P_AIR;
    _a2 = leaves.PSM.a_net;

    return (_a2 - _a1) / (_e2 - _e1)
);

∂A∂E(leaves::Leaves2D{FT}, air::AirLayer{FT}) where {FT<:AbstractFloat} = (
    (; P_AIR) = air;

    # compute the A and E at the current setting
    _gs1 = leaves.g_H₂O_s_shaded;
    _gh1 = 1 / (1 / _gs1 + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e1  = _gh1 * (saturation_vapor_pressure(leaves.t) - air.p_H₂O) / P_AIR;
    _a1  = leaves.a_net_shaded;

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    _gs2 = _gs1 + FT(0.0001);
    _gh2 = 1 / (1 / _gs2 + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _gc2 = 1 / (FT(1.6) / _gs2 + 1 / leaves.g_CO₂_b);
    leaf_photosynthesis!(leaves, air, _gc2, leaves.ppar_shaded, leaves.t);
    _e2 = _gh2 * (saturation_vapor_pressure(leaves.t) - air.p_H₂O) / P_AIR;
    _a2 = leaves.PSM.a_net;

    return (_a2 - _a1) / (_e2 - _e1)
);

∂A∂E(leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int) where {FT<:AbstractFloat} = (
    (; P_AIR) = air;

    # compute the A and E at the current setting
    _gs1 = leaves.g_H₂O_s_sunlit[ind];
    _gh1 = 1 / (1 / _gs1 + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e1  = _gh1 * (saturation_vapor_pressure(leaves.t) - air.p_H₂O) / P_AIR;
    _a1  = leaves.a_net_sunlit[ind];

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    _gs2 = _gs1 + FT(0.0001);
    _gh2 = 1 / (1 / _gs2 + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _gc2 = 1 / (FT(1.6) / _gs2 + 1 / leaves.g_CO₂_b);
    leaf_photosynthesis!(leaves, air, _gc2, leaves.ppar_sunlit[ind], leaves.t);
    _e2 = _gh2 * (saturation_vapor_pressure(leaves.t) - air.p_H₂O) / P_AIR;
    _a2 = leaves.PSM.a_net;

    return (_a2 - _a1) / (_e2 - _e1)
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-11: add function to compute ∂R∂E
#
#######################################################################################################################################################################################################
"""

    ∂R∂E(lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}) where {FT<:AbstractFloat}

Returns the marginal increase in leaf respiration rate per transpiration rate, given
- `lf` `Leaf`, `Leaves1D`, or `Leaves2D` type leaf
- `air` `AirLayer` type environmental conditions

"""
function ∂R∂E end

∂R∂E(lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}) where {FT<:AbstractFloat} = ∂R∂E(lf.SM, lf, air);

∂R∂E(sm::WangSM{FT}, lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}) where {FT<:AbstractFloat} = ∂R∂T(lf) * ∂T∂E(lf, air, sm.f_view);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-11: add function to compute ∂T∂E
#
#######################################################################################################################################################################################################
"""

    ∂T∂E(lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}, f_view::FT) where {FT<:AbstractFloat}

Returns the marginal increase in leaf temperature per transpiration rate, given
- `lf` `Leaf`, `Leaves1D`, or `Leaves2D` type leaf
- `air` `AirLayer` type environmental conditions
- `f_view` Ratio that leaf area is exposed to external sources/sinks (not other leaves, e.g., 2/LAI for canopy on average)

"""
function ∂T∂E end

∂T∂E(lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}, f_view::FT) where {FT<:AbstractFloat} = ∂T∂E(lf.BIO, lf, air, f_view);

∂T∂E(bio::BroadbandLeafBiophysics{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, f_view::FT) where {FT<:AbstractFloat} = ∂T∂E(f_view, leaf.t, leaf.WIDTH, air.wind, bio.ϵ_LW);

∂T∂E(bio::HyperspectralLeafBiophysics{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, f_view::FT) where {FT<:AbstractFloat} = ∂T∂E(f_view, leaf.t, leaf.WIDTH, air.wind, 1 - bio.τ_LW);

∂T∂E(bio::BroadbandLeafBiophysics{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, f_view::FT) where {FT<:AbstractFloat} = ∂T∂E(f_view, leaves.t[1], leaves.WIDTH, air.wind, bio.ϵ_LW);

∂T∂E(bio::HyperspectralLeafBiophysics{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, f_view::FT) where {FT<:AbstractFloat} = ∂T∂E(f_view, leaves.t[1], leaves.WIDTH, air.wind, 1 - bio.τ_LW);

∂T∂E(bio::BroadbandLeafBiophysics{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, f_view::FT) where {FT<:AbstractFloat} = ∂T∂E(f_view, leaves.t, leaves.WIDTH, air.wind, bio.ϵ_LW);

∂T∂E(bio::HyperspectralLeafBiophysics{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, f_view::FT) where {FT<:AbstractFloat} = ∂T∂E(f_view, leaves.t, leaves.WIDTH, air.wind, 1 - bio.τ_LW);

∂T∂E(f_view::FT, t::FT, width::FT, wind::FT, ϵ::FT) where {FT<:AbstractFloat} = (
    _λ = latent_heat_vapor(t) * M_H₂O(FT);
    _g = FT(1.4) * FT(0.135) * sqrt(wind / (FT(0.72) * width));
    _d = 2 * CP_D_MOL(FT) * _g + 4 * f_view * K_STEFAN(FT) * ϵ * t ^ 3;

    return _λ / _d
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-07: migrate function from older version
#     2022-Jul-11: add docs
#
#######################################################################################################################################################################################################
"""
This function returns the marginal risk for stomatal opening. This function supports a variety of optimality models for
- Leaf
- Leaves1D (ind=1 for sunlit leaves, ind=2 for shaded leaves)
- Leaves2D (ind=NA for shaded leaves, ind>1 for sunlit leaves)

"""
function ∂Θ∂E end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-08: add method for AndereggSM model on Leaf
#     2022-Jul-08: add method for EllerSM model on Leaf
#     2022-Jul-08: add method for SperrySM model on Leaf
#     2022-Jul-08: add method for WangSM model on Leaf
#     2022-Jul-08: add method for Wang2SM model on Leaf
#     2022-Jul-08: use ∂E∂P from PlantHydraulics.jl
#
#######################################################################################################################################################################################################
"""

    ∂Θ∂E(sm::AndereggSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
    ∂Θ∂E(sm::EllerSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
    ∂Θ∂E(sm::SperrySM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
    ∂Θ∂E(sm::WangSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
    ∂Θ∂E(sm::Wang2SM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}

Return the marginal risk for stomatal opening, given
- `sm` `AndereggSM`, `EllerSM`, `SperrySM`, `WangSM`, or `Wang2SM` type optimality model
- `leaf` `Leaf` type struct
- `air` `AirLayer` for environmental conditions
- `δe` Incremental flow rate to compute ∂E∂P

"""
∂Θ∂E(sm::AndereggSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; A, B) = sm;
    (; HS) = leaf;
    (; P_AIR) = air;

    # compute the E at the current setting
    _gs = leaf.g_H₂O_s;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaf.g_CO₂_b));
    _e  = _gh * (saturation_vapor_pressure(leaf.t) - air.p_H₂O) / P_AIR;

    _∂E∂P = ∂E∂P(leaf, _e; δe = δe);

    return (-2 * A * HS._p_element[end] + B) / _∂E∂P
);

∂Θ∂E(sm::EllerSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; HS) = leaf;
    (; P_AIR) = air;

    # compute the E at the current setting
    _gs = leaf.g_H₂O_s;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaf.g_CO₂_b));
    _e  = _gh * (saturation_vapor_pressure(leaf.t) - air.p_H₂O) / P_AIR;

    _∂E∂P_1 = ∂E∂P(leaf, _e; δe = δe);
    _∂E∂P_2 = ∂E∂P(leaf, _e; δe = -δe);
    _∂K∂E   = (_∂E∂P_2 - _∂E∂P_1) / δe;

    return _∂K∂E * leaf.a_net / _∂E∂P_1
);

∂Θ∂E(sm::SperrySM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; HS) = leaf;
    (; P_AIR) = air;

    # compute the E at the current setting
    _gs = leaf.g_H₂O_s;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaf.g_CO₂_b));
    _e  = _gh * (saturation_vapor_pressure(leaf.t) - air.p_H₂O) / P_AIR;

    _∂E∂P_1 = ∂E∂P(leaf, _e; δe = δe);
    _∂E∂P_2 = ∂E∂P(leaf, _e; δe = -δe);
    _∂E∂P_m = ∂E∂P(leaf, FT(0); δe = δe);
    _∂K∂E   = (_∂E∂P_2 - _∂E∂P_1) / δe;

    # compute maximum A
    _ghm = HS._e_crit / (saturation_vapor_pressure(leaf.t) - air.p_H₂O) * P_AIR;
    _gsm = 1 / (1 / _ghm - 1 / (FT(1.35) * leaf.g_CO₂_b));
    _gcm = 1 / (FT(1.6) / _gsm + 1 / leaf.g_CO₂_b);
    leaf_photosynthesis!(leaf, air, _gcm, leaf.ppar, leaf.t);
    _am = leaf.PSM.a_net;

    return _∂K∂E * _am / _∂E∂P_m
);

∂Θ∂E(sm::WangSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; HS) = leaf;
    (; P_AIR) = air;

    # compute the A and E at the current setting
    _gs = leaf.g_H₂O_s;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaf.g_CO₂_b));
    _e  = _gh * (saturation_vapor_pressure(leaf.t) - air.p_H₂O) / P_AIR;

    return leaf.a_net / max(eps(FT), (HS._e_crit - _e))
);

∂Θ∂E(sm::Wang2SM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; A) = sm;
    (; HS) = leaf;
    (; P_AIR) = air;

    # compute the E at the current setting
    _gs = leaf.g_H₂O_s;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaf.g_CO₂_b));
    _e  = _gh * (saturation_vapor_pressure(leaf.t) - air.p_H₂O) / P_AIR;

    _∂E∂P = ∂E∂P(leaf, _e; δe = δe);

    return (-1 * A * HS._p_element[end] * leaf.a_net) / _∂E∂P
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-11: add method for AndereggSM model on Leaves1D
#     2022-Jul-11: add method for EllerSM model on Leaves1D
#     2022-Jul-11: add method for SperrySM model on Leaves1D
#     2022-Jul-11: add method for WangSM model on Leaves1D
#     2022-Jul-11: add method for Wang2SM model on Leaves1D
#
#######################################################################################################################################################################################################
"""

    ∂Θ∂E(sm::AndereggSM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
    ∂Θ∂E(sm::EllerSM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
    ∂Θ∂E(sm::SperrySM{FT}, leaf::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
    ∂Θ∂E(sm::WangSM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
    ∂Θ∂E(sm::Wang2SM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}

Return the marginal risk for stomatal opening, given
- `sm` `AndereggSM`, `EllerSM`, `SperrySM`, `WangSM`, or `Wang2SM` type optimality model
- `leaves` `Leaves1D` type struct
- `air` `AirLayer` for environmental conditions
- `ind` Leaf index (1 for sunlit and 2 for shaded)
- `δe` Incremental flow rate to compute ∂E∂P

"""
∂Θ∂E(sm::AndereggSM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; A, B) = sm;
    (; HS, HS2) = leaves;
    (; P_AIR) = air;

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b[ind]));
    _e  = _gh * (saturation_vapor_pressure(leaves.t[ind]) - air.p_H₂O) / P_AIR;

    _∂E∂P = ∂E∂P(leaves, _e, ind; δe = δe);

    _hs = (ind == 1 ? HS : HS2);

    return (-2 * A * _hs._p_element[end] + B) / _∂E∂P
);

∂Θ∂E(sm::EllerSM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; HS, HS2) = leaves;
    (; P_AIR) = air;

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b[ind]));
    _e  = _gh * (saturation_vapor_pressure(leaves.t[ind]) - air.p_H₂O) / P_AIR;

    _∂E∂P_1 = ∂E∂P(leaves, _e, ind; δe = δe);
    _∂E∂P_2 = ∂E∂P(leaves, _e, ind; δe = -δe);
    _∂K∂E   = (_∂E∂P_2 - _∂E∂P_1) / δe;

    return _∂K∂E * leaves.a_net[ind] / _∂E∂P_1
);

∂Θ∂E(sm::SperrySM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; HS, HS2) = leaves;
    (; P_AIR) = air;

    _hs = (ind == 1 ? HS : HS2);

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b[ind]));
    _e  = _gh * (saturation_vapor_pressure(leaves.t[ind]) - air.p_H₂O) / P_AIR;

    _∂E∂P_1 = ∂E∂P(leaves, _e, ind; δe = δe);
    _∂E∂P_2 = ∂E∂P(leaves, _e, ind; δe = -δe);
    _∂E∂P_m = ∂E∂P(leaves, FT(0), ind; δe = δe);
    _∂K∂E   = (_∂E∂P_2 - _∂E∂P_1) / δe;

    # compute maximum A
    _ghm = _hs._e_crit / (saturation_vapor_pressure(leaves.t[ind]) - air.p_H₂O) * P_AIR;
    _gsm = 1 / (1 / _ghm - 1 / (FT(1.35) * leaves.g_CO₂_b[ind]));
    _gcm = 1 / (FT(1.6) / _gsm + 1 / leaves.g_CO₂_b[ind]);
    leaf_photosynthesis!(leaves, air, _gcm, leaves.ppar[ind], leaves.t[ind]);
    _am = leaves.PSM.a_net;

    return _∂K∂E * _am / _∂E∂P_m
);

∂Θ∂E(sm::WangSM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; HS, HS2) = leaves;
    (; P_AIR) = air;

    _hs = (ind == 1 ? HS : HS2);

    # compute the A and E at the current setting
    _gs = leaves.g_H₂O_s[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b[ind]));
    _e  = _gh * (saturation_vapor_pressure(leaves.t[ind]) - air.p_H₂O) / P_AIR;

    return leaves.a_net[ind] / max(eps(FT), (HS._e_crit - _e))
);

∂Θ∂E(sm::Wang2SM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; A) = sm;
    (; HS, HS2) = leaves;
    (; P_AIR) = air;

    _hs = (ind == 1 ? HS : HS2);

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b[ind]));
    _e  = _gh * (saturation_vapor_pressure(leaves.t[ind]) - air.p_H₂O) / P_AIR;

    _∂E∂P = ∂E∂P(leaves, _e, ind; δe = δe);

    return (-1 * A * _hs._p_element[end] * leaves.a_net[ind]) / _∂E∂P
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-11: add method for AndereggSM model on Leaves2D for shaded leaves
#     2022-Jul-11: add method for EllerSM model on Leaves2D for shaded leaves
#     2022-Jul-11: add method for SperrySM model on Leaves2D for shaded leaves
#     2022-Jul-11: add method for WangSM model on Leaves2D for shaded leaves
#     2022-Jul-11: add method for Wang2SM model on Leaves2D for shaded leaves
#
#######################################################################################################################################################################################################
"""

    ∂Θ∂E(sm::AndereggSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
    ∂Θ∂E(sm::EllerSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
    ∂Θ∂E(sm::SperrySM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
    ∂Θ∂E(sm::WangSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
    ∂Θ∂E(sm::Wang2SM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}

Return the marginal risk for stomatal opening, given
- `sm` `AndereggSM`, `EllerSM`, `SperrySM`, `WangSM`, or `Wang2SM` type optimality model
- `leaves` `Leaves2D` type struct
- `air` `AirLayer` for environmental conditions
- `δe` Incremental flow rate to compute ∂E∂P

"""
∂Θ∂E(sm::AndereggSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; A, B) = sm;
    (; HS) = leaves;
    (; P_AIR) = air;

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s_shaded;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * (saturation_vapor_pressure(leaves.t) - air.p_H₂O) / P_AIR;

    _∂E∂P = ∂E∂P(leaves, _e; δe = δe);

    return (-2 * A * HS._p_element[end] + B) / _∂E∂P
);

∂Θ∂E(sm::EllerSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; HS) = leaves;
    (; P_AIR) = air;

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s_shaded;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * (saturation_vapor_pressure(leaves.t) - air.p_H₂O) / P_AIR;

    _∂E∂P_1 = ∂E∂P(leaves, _e; δe = δe);
    _∂E∂P_2 = ∂E∂P(leaves, _e; δe = -δe);
    _∂K∂E   = (_∂E∂P_2 - _∂E∂P_1) / δe;

    return _∂K∂E * leaves.a_net_shaded / _∂E∂P_1
);

∂Θ∂E(sm::SperrySM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; HS) = leaves;
    (; P_AIR) = air;

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s_shaded;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * (saturation_vapor_pressure(leaves.t) - air.p_H₂O) / P_AIR;

    _∂E∂P_1 = ∂E∂P(leaves, _e; δe = δe);
    _∂E∂P_2 = ∂E∂P(leaves, _e; δe = -δe);
    _∂E∂P_m = ∂E∂P(leaves, FT(0); δe = δe);
    _∂K∂E   = (_∂E∂P_2 - _∂E∂P_1) / δe;

    # compute maximum A
    _ghm = HS._e_crit / (saturation_vapor_pressure(leaves.t) - air.p_H₂O) * P_AIR;
    _gsm = 1 / (1 / _ghm - 1 / (FT(1.35) * leaves.g_CO₂_b));
    _gcm = 1 / (FT(1.6) / _gsm + 1 / leaves.g_CO₂_b);
    leaf_photosynthesis!(leaves, air, _gcm, leaves.ppar_shaded, leaves.t);
    _am = leaves.PSM.a_net;

    return _∂K∂E * _am / _∂E∂P_m
);

∂Θ∂E(sm::WangSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; HS) = leaves;
    (; P_AIR) = air;

    # compute the A and E at the current setting
    _gs = leaves.g_H₂O_s_shaded;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * (saturation_vapor_pressure(leaves.t) - air.p_H₂O) / P_AIR;

    return leaves.a_net_shaded / max(eps(FT), (HS._e_crit - _e))
);

∂Θ∂E(sm::Wang2SM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; A) = sm;
    (; HS) = leaves;
    (; P_AIR) = air;

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s_shaded;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * (saturation_vapor_pressure(leaves.t) - air.p_H₂O) / P_AIR;

    _∂E∂P = ∂E∂P(leaves, _e; δe = δe);

    return (-1 * A * HS._p_element[end] * leaves.a_net_shaded) / _∂E∂P
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-11: add method for AndereggSM model on Leaves2D for sunlit leaves
#     2022-Jul-11: add method for EllerSM model on Leaves2D for sunlit leaves
#     2022-Jul-11: add method for SperrySM model on Leaves2D for sunlit leaves
#     2022-Jul-11: add method for WangSM model on Leaves2D for sunlit leaves
#     2022-Jul-11: add method for Wang2SM model on Leaves2D for sunlit leaves
#     2023-Mar-02: add a eps(FT) controller to (e_crit - e)
#
#######################################################################################################################################################################################################
"""

    ∂Θ∂E(sm::AndereggSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
    ∂Θ∂E(sm::EllerSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
    ∂Θ∂E(sm::SperrySM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
    ∂Θ∂E(sm::WangSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
    ∂Θ∂E(sm::Wang2SM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}

Return the marginal risk for stomatal opening, given
- `sm` `AndereggSM`, `EllerSM`, `SperrySM`, `WangSM`, or `Wang2SM` type optimality model
- `leaf` `Leaf` type struct
- `air` `AirLayer` for environmental conditions
- `δe` Incremental flow rate to compute ∂E∂P

"""
∂Θ∂E(sm::AndereggSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; A, B) = sm;
    (; HS) = leaves;
    (; P_AIR) = air;

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s_sunlit[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * (saturation_vapor_pressure(leaves.t) - air.p_H₂O) / P_AIR;

    _∂E∂P = ∂E∂P(leaves, _e; δe = δe);

    return (-2 * A * HS._p_element[end] + B) / _∂E∂P
);

∂Θ∂E(sm::EllerSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; HS) = leaves;
    (; P_AIR) = air;

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s_sunlit[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * (saturation_vapor_pressure(leaves.t) - air.p_H₂O) / P_AIR;

    _∂E∂P_1 = ∂E∂P(leaves, _e; δe = δe);
    _∂E∂P_2 = ∂E∂P(leaves, _e; δe = -δe);
    _∂K∂E   = (_∂E∂P_2 - _∂E∂P_1) / δe;

    return _∂K∂E * leaves.a_net_sunlit[ind] / _∂E∂P_1
);

∂Θ∂E(sm::SperrySM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; HS) = leaves;
    (; P_AIR) = air;

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s_sunlit[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * (saturation_vapor_pressure(leaves.t) - air.p_H₂O) / P_AIR;

    _∂E∂P_1 = ∂E∂P(leaves, _e; δe = δe);
    _∂E∂P_2 = ∂E∂P(leaves, _e; δe = -δe);
    _∂E∂P_m = ∂E∂P(leaves, FT(0); δe = δe);
    _∂K∂E   = (_∂E∂P_2 - _∂E∂P_1) / δe;

    # compute maximum A
    _ghm = HS._e_crit / (saturation_vapor_pressure(leaves.t) - air.p_H₂O) * P_AIR;
    _gsm = 1 / (1 / _ghm - 1 / (FT(1.35) * leaves.g_CO₂_b));
    _gcm = 1 / (FT(1.6) / _gsm + 1 / leaves.g_CO₂_b);
    leaf_photosynthesis!(leaves, air, _gcm, leaves.ppar_sunlit[ind], leaves.t);
    _am = leaves.PSM.a_net;

    return _∂K∂E * _am / _∂E∂P_m
);

∂Θ∂E(sm::WangSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; HS) = leaves;
    (; P_AIR) = air;

    # compute the A and E at the current setting
    _gs = leaves.g_H₂O_s_sunlit[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * (saturation_vapor_pressure(leaves.t) - air.p_H₂O) / P_AIR;

    return leaves.a_net_sunlit[ind] / max(eps(FT), (HS._e_crit - _e))
);

∂Θ∂E(sm::Wang2SM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    (; A) = sm;
    (; HS) = leaves;
    (; P_AIR) = air;

    # compute the E at the current setting
    _gs = leaves.g_H₂O_s_sunlit[ind];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _e  = _gh * (saturation_vapor_pressure(leaves.t) - air.p_H₂O) / P_AIR;

    _∂E∂P = ∂E∂P(leaves, _e; δe = δe);

    return (-1 * A * HS._p_element[end] * leaves.a_net_sunlit[ind]) / _∂E∂P
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-11: add function for nocturnal stomatal conductance
#
#######################################################################################################################################################################################################
"""
This function returns the ∂Θₙ∂E for nocturnal stomatal opening. Currently this function only supports WangSM which has been published for the purpose of computing nocturnal stomatal conductance.
    Supports to other optimality models will be added later when I am ready to test those.

"""
function ∂Θₙ∂E end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-11: add method for WangSM model on Leaf for nocturnal transpiration
#     2022-Jul-11: add method for WangSM model on Leaves1D for nocturnal transpiration
#     2022-Jul-11: add method for WangSM model on Leaves2D for nocturnal transpiration
#
#######################################################################################################################################################################################################
"""

    ∂Θₙ∂E(lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}) where {FT<:AbstractFloat}

Return the ∂Θ∂E for nocturnal stomatal opening, given
- `lf` `Leaf`, `Leaves1D`, or `Leaves2D` type leaf
- `air` `AirLayer` type environmental conditions

"""
∂Θₙ∂E(lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}) where {FT<:AbstractFloat} = ∂Θₙ∂E(lf.SM, lf, air);

∂Θₙ∂E(sm::WangSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}) where {FT<:AbstractFloat} = (
    (; F_FITNESS) = sm;
    (; HS) = leaf;
    (; P_AIR) = air;

    # compute the A and E at the current setting
    _gs = leaf.g_H₂O_s;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaf.g_CO₂_b));
    _gc = 1 / (FT(1.6) / _gs + 1 / leaf.g_CO₂_b);
    _e  = _gh * (saturation_vapor_pressure(leaf.t) - air.p_H₂O) / P_AIR;
    leaf_photosynthesis!(leaf, air, _gc, sm.ppar_mem, leaf.t);
    _a  = leaf.PSM.a_net;

    return _a / max(eps(FT), (HS._e_crit - _e)) * F_FITNESS
);

∂Θₙ∂E(sm::WangSM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}) where {FT<:AbstractFloat} = (
    (; F_FITNESS) = sm;
    (; HS) = leaves;
    (; P_AIR) = air;

    # compute the A and E at the current setting
    _gs = leaves.g_H₂O_s[1];
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b[1]));
    _gc = 1 / (FT(1.6) / _gs + 1 / leaves.g_CO₂_b[1]);
    _e  = _gh * (saturation_vapor_pressure(leaves.t[1]) - air.p_H₂O) / P_AIR;
    leaf_photosynthesis!(leaves, air, _gc, sm.ppar_mem, leaves.t[1]);
    _a  = leaves.PSM.a_net;

    return _a / max(eps(FT), (HS._e_crit - _e)) * F_FITNESS
);

∂Θₙ∂E(sm::WangSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}) where {FT<:AbstractFloat} = (
    (; F_FITNESS) = sm;
    (; HS) = leaves;
    (; P_AIR) = air;

    # compute the A and E at the current setting
    _gs = leaves.g_H₂O_s_shaded;
    _gh = 1 / (1 / _gs + 1 / (FT(1.35) * leaves.g_CO₂_b));
    _gc = 1 / (FT(1.6) / _gs + 1 / leaves.g_CO₂_b);
    _e  = _gh * (saturation_vapor_pressure(leaves.t) - air.p_H₂O) / P_AIR;
    leaf_photosynthesis!(leaves, air, _gc, sm.ppar_mem, leaves.t);
    _a  = leaves.PSM.a_net;

    return _a / max(eps(FT), (HS._e_crit - _e)) * F_FITNESS
);

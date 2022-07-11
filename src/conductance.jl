#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-07: add new function
#     2022-Jul-11: deflate documentations
#
#######################################################################################################################################################################################################
"""
This function returns the stomatal conductance change slope. Supported functionalities are
- Leaf
- Leaves1D (ind=1 for sunlit, ind=2 for shaded leaves)
- Leaves2D (ind=NA for shaded, ind>1 for sunlit leaves)

"""
function ∂g∂t end;


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add general method for empirical models
#     2022-Jul-07: clarify that this method is only for Leaf and shaded leaves of Leaves2D
#     2022-Jul-11: deflate documentations
#     2022-Jul-11: make this method specific for Leaf
#
#######################################################################################################################################################################################################
"""

    ∂g∂t(leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat}

Return the marginal increase of stomatal conductance, given
- `leaf` `Leaf` type struct
- `air` `AirLayer` type environmental conditions
- `β` Tuning factor (only used for empirical models)
- `δe` Incremental flow rate to compute ∂E∂P (only used for optimality models)
"""
∂g∂t(leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = ∂g∂t(leaf.SM, leaf, air; β = β, δe = δe);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    return ∂g∂t(sm, leaf, air, sm.Β.PARAM_Y; β = β)
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}, βt::BetaParameterG1; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaf, air; β = β);

    return (_gsw - leaf.g_H₂O_s) / sm.Τ
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}, βt::BetaParameterVcmax; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaf, air; β = FT(1));

    return (_gsw - leaf.g_H₂O_s) / sm.Τ
);

∂g∂t(sm::Union{AndereggSM{FT}, EllerSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    return sm.K * (∂A∂E(leaf, air) - ∂Θ∂E(sm, leaf, air; δe = δe))
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-07: add new function
#     2022-Jul-11: deflate documentations
#
#######################################################################################################################################################################################################
"""
This function updates stomatal conductance for H₂O and CO₂. Supported functionalities are
- Update conductance for H₂O prognostically
- Update conductance for CO₂ based on that for H₂O

"""
function stomatal_conductance! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add new method to update stomatal conductance prognostically
#
#######################################################################################################################################################################################################
"""

    stomatal_conductance!(leaf::Leaf{FT}, air::AirLayer{FT}, Δt::FT; β::FT = FT(1)) where {FT<:AbstractFloat}
    stomatal_conductance!(leaves::Leaves1D{FT}, air::AirLayer{FT}, Δt::FT; β::FT = FT(1)) where {FT<:AbstractFloat}
    stomatal_conductance!(leaves::Leaves2D{FT}, air::AirLayer{FT}, Δt::FT; β::FT = FT(1)) where {FT<:AbstractFloat}

Update stomatal conductance for H₂O, given
- `leaf` or `leaves` `Leaf`, `Leaves1D`, or `Leaves2D` type leaf
- `air` `AirLayer` type environmental conditions
- `Δt` Time step length `[s]`
- `β` Tuning factor for G1 or Vcmax (only for empirical models)
"""
stomatal_conductance!(leaf::Leaf{FT}, air::AirLayer{FT}, Δt::FT; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    leaf.g_H₂O_s += ∂g∂t(leaf, air; β = β) * Δt;

    return nothing
);

stomatal_conductance!(leaves::Leaves1D{FT}, air::AirLayer{FT}, Δt::FT; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    leaves.g_H₂O_s[1] += ∂g∂t(leaves, air, 1; β = β) * Δt;
    leaves.g_H₂O_s[2] += ∂g∂t(leaves, air, 2; β = β) * Δt;

    return nothing
);

stomatal_conductance!(leaves::Leaves2D{FT}, air::AirLayer{FT}, Δt::FT; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    leaves.g_H₂O_s_shaded += ∂g∂t(leaves, air; β = β) * Δt;
    for _i in eachindex(leaves.g_H₂O_s_sunlit)
        leaves.g_H₂O_s_sunlit[_i] += ∂g∂t(leaves, air, _i; β = β) * Δt;
    end;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add new method to update stomatal conductance for CO₂ based on that of H₂O
#
#######################################################################################################################################################################################################
"""

    stomatal_conductance!(leaf::Leaf{FT}) where {FT<:AbstractFloat}
    stomatal_conductance!(leaves::Leaves1D{FT}) where {FT<:AbstractFloat}
    stomatal_conductance!(leaves::Leaves2D{FT}) where {FT<:AbstractFloat}

Update stomatal conductance for CO₂ based on that of H₂O, given
- `leaf` or `leaves` `Leaf`, `Leaves1D`, or `Leaves2D` type leaf

"""
stomatal_conductance!(leaf::Leaf{FT}) where {FT<:AbstractFloat} = (
    limit_stomatal_conductance!(leaf);

    leaf.g_CO₂ = 1 / (1 / leaf.g_CO₂_b + FT(1.6) / leaf.g_H₂O_s);

    return nothing
);

stomatal_conductance!(leaves::Leaves1D{FT}) where {FT<:AbstractFloat} = (
    limit_stomatal_conductance!(leaves);

    leaves.g_CO₂[1] = 1 / (1 / leaves.g_CO₂_b[1] + FT(1.6) / leaves.g_H₂O_s[1]);
    leaves.g_CO₂[2] = 1 / (1 / leaves.g_CO₂_b[2] + FT(1.6) / leaves.g_H₂O_s[2]);

    return nothing
);

stomatal_conductance!(leaves::Leaves2D{FT}) where {FT<:AbstractFloat} = (
    limit_stomatal_conductance!(leaves);

    leaves.g_CO₂_shaded = 1 / (1 / leaves.g_CO₂_b + FT(1.6) / leaves.g_H₂O_s_shaded);
    for _i in eachindex(leaves.g_H₂O_s_sunlit)
        leaves.g_CO₂_sunlit[_i] = 1 / (1 / leaves.g_CO₂_b + FT(1.6) / leaves.g_H₂O_s_sunlit[_i]);
    end;

    return nothing
);

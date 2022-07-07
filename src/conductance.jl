#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-07: add new function
#
#######################################################################################################################################################################################################
"""
This function returns the stomatal conductance change slope. Supported methods are

$(METHODLIST)

"""
function ∂g∂t end;


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add general method for empirical models
#     2022-Jul-07: clarify that this method is only for Leaf and shaded leaves of Leaves2D
#
#######################################################################################################################################################################################################
"""

    ∂g∂t(lf::Union{Leaf{FT}, Leaves2D{FT}}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat}

A wrapper function to return the marginal increase of stomatal conductance for H₂O, given
- `leaf` `Leaf`, and `Leaves2D` (shaded fraction) type leaf
- `air` `AirLayer` type environmental conditions
- `β` Tuning factor for G1 or Vcmax
"""
∂g∂t(lf::Union{Leaf{FT}, Leaves2D{FT}}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} = ∂g∂t(lf.SM, lf, air; β = β);


∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, lf::Union{Leaf{FT}, Leaves2D{FT}}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} =
    ∂g∂t(sm, lf, air, sm.Β.PARAM_Y; β = β);


∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}, βt::BetaParameterG1; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaf, air; β = β);

    return (_gsw - leaf.g_H₂O_s) / sm.Τ
);


∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}, βt::BetaParameterVcmax; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaf, air; β = FT(1));

    return (_gsw - leaf.g_H₂O_s) / sm.Τ
);


∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, βt::BetaParameterG1; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaves, air; β = β);

    return (_gsw - leaves.g_H₂O_s_shaded) / sm.Τ
);


∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, βt::BetaParameterVcmax; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaves, air; β = FT(1));

    return (_gsw - leaves.g_H₂O_s_shaded) / sm.Τ
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add general method for empirical models
#     2022-Jul-07: clarify that this method is only for Leaves1D and sunlit leaves of Leaves2D
#
#######################################################################################################################################################################################################
"""

    ∂g∂t(lf::Union{Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat}

A wrapper function to return the marginal increase of stomatal conductance for H₂O, given
- `leaf` `Leaves1D`, and `Leaves2D` (sunlit fraction) type leaf
- `air` `AirLayer` type environmental conditions
- `ind` Index of the leaves (1 for sunlit and 2 for shaded for Leaves1D, all sunlit for Leaves2D)
- `β` Tuning factor for G1 or Vcmax
"""
∂g∂t(lf::Union{Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = ∂g∂t(lf.SM, lf, air, ind; β = β);


∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, lf::Union{Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} =
    ∂g∂t(sm, lf, air, sm.Β.PARAM_Y, ind; β = β);


∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves1D{FT}, air::AirLayer{FT}, βt::BetaParameterG1, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaves, air, ind; β = β);

    return (_gsw - leaves.g_H₂O_s[ind]) / sm.Τ
);


∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves1D{FT}, air::AirLayer{FT}, βt::BetaParameterVcmax, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaves, air, ind; β = FT(1));

    return (_gsw - leaves.g_H₂O_s[ind]) / sm.Τ
);


∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, βt::BetaParameterG1, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaves, air, ind; β = β);

    return (_gsw - leaves.g_H₂O_s_sunlit[ind]) / sm.Τ
);


∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, βt::BetaParameterVcmax, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaves, air, ind; β = FT(1));

    return (_gsw - leaves.g_H₂O_s_sunlit[ind]) / sm.Τ
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-07: add new function
#
#######################################################################################################################################################################################################
"""
This function returns the stomatal conductance change slope. Supported methods are

$(METHODLIST)

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
#     2022-Jul-07: add new method to update stomatal conductance prognostically
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

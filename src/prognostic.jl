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
#
#######################################################################################################################################################################################################
"""

    ∂g∂t(leaf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat}

A wrapper function to return the marginal increase of stomatal conductance for H₂O, given
- `leaf` `Leaf`, `Leaves1D`, and `Leaves2D` type leaf
- `air` `AirLayer` type environmental conditions
- `β` Tuning factor for G1 or Vcmax
"""
∂g∂t(leaf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} = ∂g∂t(leaf.SM, leaf, air; β = β);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add general method for empirical models
#
#######################################################################################################################################################################################################
"""

    ∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat}

A wrapper function to return the marginal increase of stomatal conductance for H₂O from empirical models, given
- `sm` `BallBerrySM`, `GentineSM`, `LeuningSM`, and `MedlynSM` type empirical models
- `leaf` `Leaf`, `Leaves1D`, and `Leaves2D` type leaf
- `air` `AirLayer` type environmental conditions
- `β` Tuning factor for G1 or Vcmax
"""
∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} =
    ∂g∂t(sm, leaf, air, sm.Β.PARAM_Y; β = β);


∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}, βt::BetaParameterG1; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaf, air; β = β);

    return (_gsw - leaf.g_H₂O_s) / sm.Τ
);


∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}, βt::BetaParameterVcmax; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaf, air; β = FT(1));

    return (_gsw - leaf.g_H₂O_s) / sm.Τ
);


∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves1D{FT}, air::AirLayer{FT}, βt::BetaParameterG1, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaves, air, ind; β = β);

    return (_gsw - leaves.g_H₂O_s[ind]) / sm.Τ
);


∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves1D{FT}, air::AirLayer{FT}, βt::BetaParameterVcmax, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaves, air, ind; β = FT(1));

    return (_gsw - leaves.g_H₂O_s[ind]) / sm.Τ
);


∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, βt::BetaParameterG1; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaves, air; β = β);

    return (_gsw - leaves.g_H₂O_s_shaded) / sm.Τ
);


∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, βt::BetaParameterVcmax; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaves, air; β = FT(1));

    return (_gsw - leaves.g_H₂O_s_shaded) / sm.Τ
);


∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, βt::BetaParameterG1, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaves, air, ind; β = β);

    return (_gsw - leaves.g_H₂O_s_sunlit[ind]) / sm.Τ
);


∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, βt::BetaParameterVcmax, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaves, air, ind; β = FT(1));

    return (_gsw - leaves.g_H₂O_s_sunlit[ind]) / sm.Τ
);



function stomatal_conductance! end



stomatal_conductance!(leaf::Leaf{FT}, air::AirLayer{FT}, Δt::FT = 60; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack SM = leaf;

    leaf.g_H₂O_s += ∂g∂t(leaf, air; β = β) * Δt;

    return nothing
);

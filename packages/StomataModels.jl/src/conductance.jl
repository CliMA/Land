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

∂g∂t(sm::Union{AndereggSM{FT}, EllerSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    return sm.K * (∂A∂E(leaf, air) - ∂Θ∂E(sm, leaf, air; δe = δe))
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    return ∂g∂t(sm, leaf, air, sm.β.PARAM_Y; β = β)
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}, βt::BetaParameterG1; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaf, air; β = β);

    return (_gsw - leaf.g_H₂O_s) / sm.τ
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}, βt::BetaParameterVcmax; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaf, air; β = FT(1));

    return (_gsw - leaf.g_H₂O_s) / sm.τ
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-11: add method for Leaves1D
#
#######################################################################################################################################################################################################
"""

    ∂g∂t(leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat}

Return the marginal increase of stomatal conductance, given
- `leaves` `Leaves1D` type struct
- `air` `AirLayer` type environmental conditions
- `ind` Leaf index (1 for sunlit and 2 for shaded)
- `β` Tuning factor (only used for empirical models)
- `δe` Incremental flow rate to compute ∂E∂P (only used for optimality models)

"""
∂g∂t(leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = ∂g∂t(leaves.SM, leaves, air, ind; β = β, δe = δe);

∂g∂t(sm::Union{AndereggSM{FT}, EllerSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    return sm.K * (∂A∂E(leaves, air, ind) - ∂Θ∂E(sm, leaves, air, ind; δe = δe))
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    return ∂g∂t(sm, leaves, air, sm.β.PARAM_Y, ind; β = β)
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves1D{FT}, air::AirLayer{FT}, βt::BetaParameterG1, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaves, air, ind; β = β);

    return (_gsw - leaves.g_H₂O_s[ind]) / sm.τ
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves1D{FT}, air::AirLayer{FT}, βt::BetaParameterVcmax, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaves, air, ind; β = FT(1));

    return (_gsw - leaves.g_H₂O_s[ind]) / sm.τ
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-11: add method for Leaves2D shaded leaves
#
#######################################################################################################################################################################################################
"""

    ∂g∂t(leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat}

Return the marginal increase of stomatal conductance, given
- `leaves` `Leaves2D` type struct
- `air` `AirLayer` type environmental conditions
- `β` Tuning factor (only used for empirical models)
- `δe` Incremental flow rate to compute ∂E∂P (only used for optimality models)

"""
∂g∂t(leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = ∂g∂t(leaves.SM, leaves, air; β = β, δe = δe);

∂g∂t(sm::Union{AndereggSM{FT}, EllerSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    return sm.K * (∂A∂E(leaves, air) - ∂Θ∂E(sm, leaves, air; δe = δe))
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    return ∂g∂t(sm, leaves, air, sm.β.PARAM_Y; β = β)
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, βt::BetaParameterG1; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaves, air; β = β);

    return (_gsw - leaves.g_H₂O_s_shaded) / sm.τ
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, βt::BetaParameterVcmax; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaves, air; β = FT(1));

    return (_gsw - leaves.g_H₂O_s_shaded) / sm.τ
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-11: add method for Leaves2D sunlit leaves
#
#######################################################################################################################################################################################################
"""

    ∂g∂t(leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat}

Return the marginal increase of stomatal conductance, given
- `leaves` `Leaves2D` type struct
- `air` `AirLayer` type environmental conditions
- `ind` Sunlit leaf index within the leaf angular distribution
- `β` Tuning factor (only used for empirical models)
- `δe` Incremental flow rate to compute ∂E∂P (only used for optimality models)

"""
∂g∂t(leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = ∂g∂t(leaves.SM, leaves, air, ind; β = β, δe = δe);

∂g∂t(sm::Union{AndereggSM{FT}, EllerSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat} =
    sm.K * (∂A∂E(leaves, air, ind) - ∂Θ∂E(sm, leaves, air, ind; δe = δe));

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    return ∂g∂t(sm, leaves, air, sm.β.PARAM_Y, ind; β = β)
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, βt::BetaParameterG1, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaves, air, ind; β = β);

    return (_gsw - leaves.g_H₂O_s_sunlit[ind]) / sm.τ
);

∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaves::Leaves2D{FT}, air::AirLayer{FT}, βt::BetaParameterVcmax, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaves, air, ind; β = FT(1));

    return (_gsw - leaves.g_H₂O_s_sunlit[ind]) / sm.τ
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-11: add method for nocturnal transpiration for WangSM model
#     2022-Jul-11: rename function to ∂gₙ∂t
#
#######################################################################################################################################################################################################
"""

    ∂gₙ∂t(lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}) where {FT<:AbstractFloat}

Return the marginal increase of stomatal conductance, given
- `lf` `Leaf`, `Leaves1D`, or `Leaves2D` type struct
- `air` `AirLayer` type environmental conditions

"""
function ∂gₙ∂t end

∂gₙ∂t(lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}) where {FT<:AbstractFloat} = ∂gₙ∂t(lf.SM, lf, air);

∂gₙ∂t(sm::WangSM{FT}, lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}) where {FT<:AbstractFloat} = sm.K * (∂R∂E(lf, air) - ∂Θₙ∂E(lf, air));


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
#     2022-Jul-12: add new method to compute marginal stomatal conductance increase
#     2022-Jul-12: add new method to update marginal stomatal conductance for SPAC
# To do
#     TODO: be careful with the β here (need to used the value stored in empirical SM)
#     TODO: use ∂gₙ∂t for nighttime conditions
#
#######################################################################################################################################################################################################
"""

    stomatal_conductance!(spac::MonoElementSPAC{FT}; β::FT = FT(1)) where {FT<:AbstractFloat}
    stomatal_conductance!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}; β::FT = FT(1)) where {FT<:AbstractFloat}

Update marginal stomatal conductance, given
- `spac` `MonoElementSPAC`, `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` type struct
- `β` Tuning factor

"""
stomatal_conductance!(spac::MonoElementSPAC{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    (; AIR, LEAF) = spac;

    stomatal_conductance!(LEAF, AIR; β = β);

    return nothing
);

stomatal_conductance!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    (; AIR, LEAVES, LEAVES_INDEX) = spac;

    for _i in eachindex(LEAVES_INDEX)
        stomatal_conductance!(LEAVES[_i], AIR[LEAVES_INDEX[_i]]; β = β);
    end;

    return nothing
);

stomatal_conductance!(leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    leaf.∂g∂t = ∂g∂t(leaf, air; β = β);

    return nothing
);

stomatal_conductance!(leaves::Leaves1D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    leaves.∂g∂t[1] = ∂g∂t(leaves, air, 1; β = β);
    leaves.∂g∂t[2] = ∂g∂t(leaves, air, 2; β = β);

    return nothing
);

stomatal_conductance!(leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    leaves.∂g∂t_shaded = ∂g∂t(leaves, air; β = β);
    for _i in eachindex(leaves.∂g∂t_sunlit)
        leaves.∂g∂t_sunlit[_i] = ∂g∂t(leaves, air, _i; β = β);
    end;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add new method to update stomatal conductance for CO₂ based on that of H₂O
#     2022-Jul-07: add new method to update stomatal conductance prognostically
#     2022-Jul-12: move ∂g∂t to another method
#     2022-Jul-12: add method to update g for SPAC
#     2022-Jul-26: limit g in range after updating stomatal conductance
#
#######################################################################################################################################################################################################
"""

    stomatal_conductance!(spac::MonoElementSPAC{FT}, Δt::FT) where {FT<:AbstractFloat}
    stomatal_conductance!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, Δt::FT) where {FT<:AbstractFloat}

Update marginal stomatal conductance, given
- `spac` `MonoElementSPAC`, `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` type struct
- `Δt` Time step length `[s]`

"""
stomatal_conductance!(spac::MonoElementSPAC{FT}, Δt::FT) where {FT<:AbstractFloat} = (
    (; LEAF) = spac;

    stomatal_conductance!(LEAF, Δt);

    return nothing
);

stomatal_conductance!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, Δt::FT) where {FT<:AbstractFloat} = (
    (; LEAVES) = spac;

    for _leaves in LEAVES
        stomatal_conductance!(_leaves, Δt);
    end;

    return nothing
);

stomatal_conductance!(leaf::Leaf{FT}, Δt::FT) where {FT<:AbstractFloat} = (
    leaf.g_H₂O_s += leaf.∂g∂t * Δt;
    stomatal_conductance!(leaf);

    return nothing
);

stomatal_conductance!(leaves::Leaves1D{FT}, Δt::FT) where {FT<:AbstractFloat} = (
    leaves.g_H₂O_s[1] += leaves.∂g∂t[1] * Δt;
    leaves.g_H₂O_s[2] += leaves.∂g∂t[2] * Δt;
    stomatal_conductance!(leaves);

    return nothing
);

stomatal_conductance!(leaves::Leaves2D{FT}, Δt::FT) where {FT<:AbstractFloat} = (
    leaves.g_H₂O_s_shaded += leaves.∂g∂t_shaded * Δt;
    for _i in eachindex(leaves.g_H₂O_s_sunlit)
        leaves.g_H₂O_s_sunlit[_i] += leaves.∂g∂t_sunlit[_i] * Δt;
    end;
    stomatal_conductance!(leaves);

    return nothing
);

stomatal_conductance!(leaf::Leaf{FT}) where {FT<:AbstractFloat} = (
    limit_stomatal_conductance!(leaf);

    leaf._g_CO₂ = 1 / (1 / leaf.g_CO₂_b + FT(1.6) / leaf.g_H₂O_s);

    return nothing
);

stomatal_conductance!(leaves::Leaves1D{FT}) where {FT<:AbstractFloat} = (
    limit_stomatal_conductance!(leaves);

    leaves._g_CO₂[1] = 1 / (1 / leaves.g_CO₂_b[1] + FT(1.6) / leaves.g_H₂O_s[1]);
    leaves._g_CO₂[2] = 1 / (1 / leaves.g_CO₂_b[2] + FT(1.6) / leaves.g_H₂O_s[2]);

    return nothing
);

stomatal_conductance!(leaves::Leaves2D{FT}) where {FT<:AbstractFloat} = (
    limit_stomatal_conductance!(leaves);

    leaves._g_CO₂_shaded = 1 / (1 / leaves.g_CO₂_b + FT(1.6) / leaves.g_H₂O_s_shaded);
    for _i in eachindex(leaves.g_H₂O_s_sunlit)
        leaves._g_CO₂_sunlit[_i] = 1 / (1 / leaves.g_CO₂_b + FT(1.6) / leaves.g_H₂O_s_sunlit[_i]);
    end;

    return nothing
);

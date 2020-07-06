###############################################################################
#
# Calculate empirical gsw from the equations
#
###############################################################################
"""
    empirical_gsw_from_model(
            model::EmpiricalStomatalModel,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            β::FT)
    empirical_gsw_from_model(
            model::EmpiricalStomatalModel,
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            β::FT)
    empirical_gsw_from_model(
            model::EmpiricalStomatalModel,
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            β::FT,
            ind::Int)

Steady state gsw from empirical approach given
- `model` [`EmpiricalStomatalModel`](@ref) type empirical model parameter set
- `leaf` [`Leaf`](@ref) type struct
- `leaves` [`Leaves`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
- `β` Correction factor over the g1 part of an empirical model
- `ind` Nth leaf in Leaves
"""
function empirical_gsw_from_model(
            model::ESMBallBerry{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            β::FT
            ) where {FT<:AbstractFloat}
    @unpack g0, g1    = model;
    @unpack An, p_s   = leaf;
    @unpack p_atm, RH = envir;

    return g0 .+ g1 * RH * p_atm * FT(1e-6) * β .* An ./ p_s
end

function empirical_gsw_from_model(
            model::ESMGentine{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            β::FT
            ) where {FT<:AbstractFloat}
    @unpack g0, g1  = model;
    @unpack An, p_i = leaf;
    @unpack p_atm   = envir;

    return g0 .+ g1 * p_atm * FT(1e-6) * β .* An ./ p_i
end

function empirical_gsw_from_model(
            model::ESMLeuning{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            β::FT
            ) where {FT<:AbstractFloat}
    @unpack d0, g0, g1             = model;
    @unpack An, p_s, p_sat, Γ_star = leaf;
    @unpack p_atm, p_H₂O           = envir;

    return g0 .+ g1 * p_atm * FT(1e-6) / (1 + (p_sat - p_H₂O)/d0) *
                 β .* An ./ (p_s .- Γ_star)
end

function empirical_gsw_from_model(
            model::ESMMedlyn{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            β::FT
            ) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack An, p_sat         = leaf;
    @unpack p_a, p_atm, p_H₂O = envir;

    return g0 .+ p_atm * FT(1e-6) / p_a * (1 + g1/sqrt(p_sat - p_H₂O)) *
                 β .* An
end

function empirical_gsw_from_model(
            model::ESMBallBerry{FT},
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            β::FT
            ) where {FT<:AbstractFloat}
    @unpack g0, g1    = model;
    @unpack An, p_s   = leaves;
    @unpack p_atm, RH = envir;

    return g0 .+ g1 * RH * p_atm * FT(1e-6) * β .* An ./ p_s
end

function empirical_gsw_from_model(
            model::ESMGentine{FT},
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            β::FT
            ) where {FT<:AbstractFloat}
    @unpack g0, g1  = model;
    @unpack An, p_i = leaves;
    @unpack p_atm   = envir;

    return g0 .+ g1 * p_atm * FT(1e-6) * β .* An ./ p_i
end

function empirical_gsw_from_model(
            model::ESMLeuning{FT},
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            β::FT
            ) where {FT<:AbstractFloat}
    @unpack d0, g0, g1             = model;
    @unpack An, p_s, p_sat, Γ_star = leaves;
    @unpack p_atm, p_H₂O           = envir;

    return g0 .+ g1 * p_atm * FT(1e-6) / (1 + (p_sat - p_H₂O)/d0) *
                 β .* An ./ (p_s .- Γ_star)
end

function empirical_gsw_from_model(
            model::ESMMedlyn{FT},
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            β::FT
            ) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack An, p_sat         = leaves;
    @unpack p_a, p_atm, p_H₂O = envir;

    return g0 .+ p_atm * FT(1e-6) / p_a * (1 + g1/sqrt(p_sat - p_H₂O)) *
                 β .* An
end

function empirical_gsw_from_model(
            model::ESMBallBerry{FT},
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            β::FT,
            ind::Int
            ) where {FT<:AbstractFloat}
    @unpack g0, g1    = model;
    @unpack An, p_s   = leaves;
    @unpack p_atm, RH = envir;

    return g0 + g1 * RH * p_atm * FT(1e-6) * β * An[ind] / p_s[ind]
end

function empirical_gsw_from_model(
            model::ESMGentine{FT},
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            β::FT,
            ind::Int
            ) where {FT<:AbstractFloat}
    @unpack g0, g1  = model;
    @unpack An, p_i = leaves;
    @unpack p_atm   = envir;

    return g0 + g1 * p_atm * FT(1e-6) * β * An[ind] / p_i[ind]
end

function empirical_gsw_from_model(
            model::ESMLeuning{FT},
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            β::FT,
            ind::Int
            ) where {FT<:AbstractFloat}
    @unpack d0, g0, g1             = model;
    @unpack An, p_s, p_sat, Γ_star = leaves;
    @unpack p_atm, p_H₂O           = envir;

    return g0 + g1 * p_atm * FT(1e-6) / (1 + (p_sat - p_H₂O)/d0) *
                β * An[ind] / (p_s[ind] - Γ_star)
end

function empirical_gsw_from_model(
            model::ESMMedlyn{FT},
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            β::FT,
            ind::Int
            ) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack An, p_sat         = leaves;
    @unpack p_a, p_atm, p_H₂O = envir;

    return g0 + p_atm * FT(1e-6) / p_a * (1 + g1/sqrt(p_sat - p_H₂O)) *
                β * An[ind]
end

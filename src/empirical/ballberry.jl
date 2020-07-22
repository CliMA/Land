###############################################################################
#
# Calculate empirical gsw from the equations
#
###############################################################################
"""
    empirical_gsw_from_model(model::EmpiricalStomatalModel{FT}, leaf::Leaf{FT}, envir::AirLayer{FT}, β::FT) where {FT<:AbstractFloat}
    empirical_gsw_from_model(model::EmpiricalStomatalModel{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT}, β::FT) where {FT<:AbstractFloat}
    empirical_gsw_from_model(model::EmpiricalStomatalModel{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT}, β::FT, ind::Int) where {FT<:AbstractFloat}

Steady state gsw from empirical approach given
- `model` [`EmpiricalStomatalModel`](@ref) type empirical model parameter set
- `leaf` [`Leaf`] type struct
- `canopyi` [`CanopyLayer`](@ref) type struct
- `envir` [`AirLayer`] type struct
- `β` Correction factor over the g1 part of an empirical model
- `ind` Nth leaf in the canopy layer
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

    return g0 + g1 * RH * p_atm * FT(1e-6) * β * An / p_s
end




function empirical_gsw_from_model(
            model::ESMBallBerry{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            β::FT
) where {FT<:AbstractFloat}
    @unpack g0, g1    = model;
    @unpack An, p_s   = canopyi;
    @unpack p_atm, RH = envir;

    return g0 .+ g1 * RH * p_atm * FT(1e-6) * β .* An ./ p_s
end




function empirical_gsw_from_model(
            model::ESMBallBerry{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            β::FT,
            ind::Int
) where {FT<:AbstractFloat}
    @unpack g0, g1    = model;
    @unpack An, p_s   = canopyi;
    @unpack p_atm, RH = envir;

    return g0 + g1 * RH * p_atm * FT(1e-6) * β * An[ind] / p_s[ind]
end

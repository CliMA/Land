###############################################################################
#
# Calculate empirical gsw from the equations
#
###############################################################################
"""
    stomatal_conductance(
                model::EmpiricalStomatalModel{FT},
                leaf::Leaf{FT},
                envir::AirLayer{FT},
                β::FT
    ) where {FT<:AbstractFloat}
    stomatal_conductance(
                model::EmpiricalStomatalModel{FT},
                canopyi::CanopyLayer{FT},
                envir::AirLayer{FT},
                β::FT
    ) where {FT<:AbstractFloat}
    stomatal_conductance(
                model::EmpiricalStomatalModel{FT},
                canopyi::CanopyLayer{FT},
                envir::AirLayer{FT},
                β::FT,
                ind::Int
    ) where {FT<:AbstractFloat}

Steady state gsw from empirical approach given
- `model` [`EmpiricalStomatalModel`](@ref) type empirical model parameter set
- `leaf` [`Leaf`] type struct
- `canopyi` [`CanopyLayer`](@ref) type struct
- `envir` [`AirLayer`] type struct
- `β` Correction factor over the g1 part of an empirical model
- `ind` Nth leaf in the canopy layer
"""
function stomatal_conductance(
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




function stomatal_conductance(
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




function stomatal_conductance(
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




function stomatal_conductance(
            model::ESMGentine{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            β::FT
) where {FT<:AbstractFloat}
    @unpack g0, g1  = model;
    @unpack An, p_i = leaf;
    @unpack p_atm   = envir;

    return g0 + g1 * p_atm * FT(1e-6) * β * An / p_i
end




function stomatal_conductance(
            model::ESMGentine{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            β::FT
) where {FT<:AbstractFloat}
    @unpack g0, g1  = model;
    @unpack An, p_i = canopyi;
    @unpack p_atm   = envir;

    return g0 .+ g1 * p_atm * FT(1e-6) * β .* An ./ p_i
end




function stomatal_conductance(
            model::ESMGentine{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            β::FT,
            ind::Int
) where {FT<:AbstractFloat}
    @unpack g0, g1  = model;
    @unpack An, p_i = canopyi;
    @unpack p_atm   = envir;

    return g0 + g1 * p_atm * FT(1e-6) * β * An[ind] / p_i[ind]
end




function stomatal_conductance(
            model::ESMLeuning{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            β::FT
) where {FT<:AbstractFloat}
    @unpack d0, g0, g1             = model;
    @unpack An, p_s, p_sat, Γ_star = leaf;
    @unpack p_atm, p_H₂O           = envir;

    return g0 + g1 * p_atm * FT(1e-6) / (1 + (p_sat - p_H₂O)/d0) *
                β * An / (p_s - Γ_star)
end




function stomatal_conductance(
            model::ESMLeuning{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            β::FT
) where {FT<:AbstractFloat}
    @unpack d0, g0, g1     = model;
    @unpack An, p_s, p_sat = canopyi;
    @unpack Γ_star         = canopyi.ps;
    @unpack p_atm, p_H₂O   = envir;

    return g0 .+ g1 * p_atm * FT(1e-6) / (1 + (p_sat - p_H₂O)/d0) *
                 β .* An ./ (p_s .- Γ_star)
end




function stomatal_conductance(
            model::ESMLeuning{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            β::FT,
            ind::Int
) where {FT<:AbstractFloat}
    @unpack d0, g0, g1     = model;
    @unpack An, p_s, p_sat = canopyi;
    @unpack Γ_star         = canopyi.ps;
    @unpack p_atm, p_H₂O   = envir;

    return g0 + g1 * p_atm * FT(1e-6) / (1 + (p_sat - p_H₂O)/d0) *
                β * An[ind] / (p_s[ind] - Γ_star)
end




function stomatal_conductance(
            model::ESMMedlyn{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            β::FT
) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack An, p_sat         = leaf;
    @unpack p_a, p_atm, p_H₂O = envir;
    vpd = max(FT(0.001), p_sat - p_H₂O);

    return g0 + p_atm * FT(1e-6) / p_a * (1 + g1/sqrt(vpd)) * β * An
end




function stomatal_conductance(
            model::ESMMedlyn{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            β::FT
) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack An, p_sat         = canopyi;
    @unpack p_a, p_atm, p_H₂O = envir;
    vpd = max(FT(0.001), p_sat - p_H₂O);

    return g0 .+ p_atm * FT(1e-6) / p_a * (1 + g1/sqrt(vpd)) * β .* An
end




function stomatal_conductance(
            model::ESMMedlyn{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            β::FT,
            ind::Int
) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack An, p_sat         = canopyi;
    @unpack p_a, p_atm, p_H₂O = envir;
    vpd = max(FT(0.001), p_sat - p_H₂O);

    return g0 + p_atm * FT(1e-6) / p_a * (1 + g1/sqrt(vpd)) * β * An[ind]
end

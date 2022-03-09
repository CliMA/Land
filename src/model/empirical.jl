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
    @unpack P_AIR, rh = envir;

    return g0 + g1 * rh * P_AIR * FT(1e-6) * β * leaf.PSM.a_net / leaf.p_CO₂_s
end




function stomatal_conductance(
            model::ESMBallBerry{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            β::FT
) where {FT<:AbstractFloat}
    @unpack g0, g1    = model;
    @unpack P_AIR, rh = envir;

    return g0 .+ g1 * rh * P_AIR * FT(1e-6) * β .* canopyi.ps.PSM.a_net / canopyi.ps.p_CO₂_s
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
    @unpack P_AIR, rh = envir;

    return g0 + g1 * rh * P_AIR * FT(1e-6) * β * An[ind] / p_s[ind]
end




function stomatal_conductance(
            model::ESMGentine{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            β::FT
) where {FT<:AbstractFloat}
    @unpack g0, g1  = model;
    @unpack P_AIR   = envir;

    return g0 + g1 * P_AIR * FT(1e-6) * β * leaf.PSM.a_net / leaf.p_CO₂_i
end




function stomatal_conductance(
            model::ESMGentine{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            β::FT
) where {FT<:AbstractFloat}
    @unpack g0, g1  = model;
    @unpack P_AIR   = envir;

    return g0 .+ g1 * P_AIR * FT(1e-6) * β .* canopyi.ps.PSM.a_net / canopyi.ps.p_CO₂_i
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
    @unpack P_AIR   = envir;

    return g0 + g1 * P_AIR * FT(1e-6) * β * An[ind] / p_i[ind]
end




function stomatal_conductance(
            model::ESMLeuning{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            β::FT
) where {FT<:AbstractFloat}
    @unpack d0, g0, g1 = model;
    @unpack P_AIR, p_H₂O = envir;

    if typeof(leaf.PSM) <: C4VJPModel
        _γ_s = FT(0);
    else
        _γ_s = leaf.PSM.γ_star;
    end;

    return g0 + g1 * P_AIR * FT(1e-6) / (1 + (leaf.p_H₂O_sat - p_H₂O)/d0) * β * leaf.PSM.a_net / (leaf.p_CO₂_s - _γ_s)
end




function stomatal_conductance(
            model::ESMLeuning{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            β::FT
) where {FT<:AbstractFloat}
    @unpack d0, g0, g1     = model;
    @unpack P_AIR, p_H₂O   = envir;

    if typeof(canopyi.ps.PSM) <: C4VJPModel
        _γ_s = FT(0);
    else
        _γ_s = canopyi.ps.PSM.γ_star;
    end;

    return g0 .+ g1 * P_AIR * FT(1e-6) / (1 + (canopyi.ps.p_H₂O_sat - p_H₂O)/d0) * β .* canopyi.ps.PSM.a_net ./ (canopyi.ps.p_CO₂_s .- _γ_s)
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
    @unpack P_AIR, p_H₂O   = envir;

    if typeof(canopyi.ps.PSM) <: C4VJPModel
        _γ_s = FT(0);
    else
        _γ_s = canopyi.ps.PSM.γ_star;
    end;

    return g0 + g1 * P_AIR * FT(1e-6) / (1 + (p_sat - p_H₂O)/d0) * β * An[ind] / (p_s[ind] - _γ_s)
end




function stomatal_conductance(
            model::ESMMedlyn{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            β::FT
) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack p_CO₂, P_AIR, p_H₂O = envir;
    vpd = max(FT(0.001), leaf.p_H₂O_sat - p_H₂O);

    return g0 + P_AIR * FT(1e-6) / p_CO₂ * (1 + g1/sqrt(vpd)) * β * leaf.PSM.a_net
end




function stomatal_conductance(
            model::ESMMedlyn{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            β::FT
) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack p_CO₂, P_AIR, p_H₂O = envir;
    vpd = max(FT(0.001), canopyi.ps.p_H₂O_sat - p_H₂O);

    return g0 .+ P_AIR * FT(1e-6) / p_CO₂ * (1 + g1/sqrt(vpd)) * β .* canopyi.ps.PSM.a_net
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
    @unpack p_CO₂, P_AIR, p_H₂O = envir;
    vpd = max(FT(0.001), p_sat - p_H₂O);

    return g0 + P_AIR * FT(1e-6) / p_CO₂ * (1 + g1/sqrt(vpd)) * β * An[ind]
end

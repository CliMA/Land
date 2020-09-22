###############################################################################
#
# Leuning model
#
###############################################################################
function empirical_gsw_from_model(
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




function empirical_gsw_from_model(
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




function empirical_gsw_from_model(
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

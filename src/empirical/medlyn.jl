###############################################################################
#
# Medlyn model
#
###############################################################################
function empirical_gsw_from_model(
            model::ESMMedlyn{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            β::FT
) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack An, p_sat         = leaf;
    @unpack p_a, p_atm, p_H₂O = envir;

    return g0 + p_atm * FT(1e-6) / p_a * (1 + g1/sqrt(p_sat - p_H₂O)) *
                β * An
end




function empirical_gsw_from_model(
            model::ESMMedlyn{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            β::FT
) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack An, p_sat         = canopyi;
    @unpack p_a, p_atm, p_H₂O = envir;

    return g0 .+ p_atm * FT(1e-6) / p_a * (1 + g1/sqrt(p_sat - p_H₂O)) *
                 β .* An
end




function empirical_gsw_from_model(
            model::ESMMedlyn{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            β::FT,
            ind::Int
) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack An, p_sat         = canopyi;
    @unpack p_a, p_atm, p_H₂O = envir;

    return g0 + p_atm * FT(1e-6) / p_a * (1 + g1/sqrt(p_sat - p_H₂O)) *
                β * An[ind]
end

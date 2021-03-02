###############################################################################
#
# Gentine model
#
###############################################################################
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

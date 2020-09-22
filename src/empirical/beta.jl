###############################################################################
#
# Beta functions to make soil moisture correction
#
###############################################################################
"""
    β_factor(bt::AbstractBetaFunction{FT},
             pl::FT,
             ps::FT,
             swc::FT) where {FT<:AbstractFloat}

Calculate the β correction factor, given
- `bt` [`AbstractBetaFunction`](@ref) type struct
- `pl` Leaf water potential `[MPa]`
- `ps` Soil water potential `[MPa]`
- `swc` Soil water content
"""
function β_factor(
            bt::BetaGLinearPleaf{FT},
            pl::FT,
            ps::FT,
            swc::FT
) where {FT<:AbstractFloat}
    @unpack p_max, p_min = bt;

    if pl >= p_max
        return FT(1)
    elseif pl <= p_min
        return FT(0)
    else
        return (pl - p_min) / (p_max - p_min)
    end
end




function β_factor(
            bt::BetaGLinearPsoil{FT},
            pl::FT,
            ps::FT,
            swc::FT
) where {FT<:AbstractFloat}
    @unpack p_max, p_min = bt;

    if ps >= p_max
        return FT(1)
    elseif ps <= p_min
        return FT(0)
    else
        return (ps - p_min) / (p_max - p_min)
    end
end




function β_factor(
            bt::BetaGLinearSWC{FT},
            pl::FT,
            ps::FT,
            swc::FT
) where {FT<:AbstractFloat}
    @unpack swc_max, swc_min = bt;

    if swc >= swc_max
        return FT(1)
    elseif swc <= swc_min
        return FT(0)
    else
        return (swc - swc_min) / (swc_max - swc_min)
    end
end




function β_factor(
            bt::BetaVLinearPleaf{FT},
            pl::FT,
            ps::FT,
            swc::FT
) where {FT<:AbstractFloat}
    @unpack p_max, p_min = bt;

    if pl >= p_max
        return FT(1)
    elseif pl <= p_min
        return FT(0)
    else
        return (pl - p_min) / (p_max - p_min)
    end
end




function β_factor(
            bt::BetaVLinearPsoil{FT},
            pl::FT,
            ps::FT,
            swc::FT
) where {FT<:AbstractFloat}
    @unpack p_max, p_min = bt;

    if ps >= p_max
        return FT(1)
    elseif ps <= p_min
        return FT(0)
    else
        return (ps - p_min) / (p_max - p_min)
    end
end




function β_factor(
            bt::BetaVLinearSWC{FT},
            pl::FT,
            ps::FT,
            swc::FT
) where {FT<:AbstractFloat}
    @unpack swc_max, swc_min = bt;

    if swc >= swc_max
        return FT(1)
    elseif swc <= swc_min
        return FT(0)
    else
        return (swc - swc_min) / (swc_max - swc_min)
    end
end

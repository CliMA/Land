###############################################################################
#
# Beta functions to make soil moisture correction
#
###############################################################################
"""
    β_factor(bt::AbstractBetaFunction{FT},
             p_leaf::FT,
             p_soil::FT,
             swc::FT
    ) where {FT<:AbstractFloat}

Calculate the β correction factor, given
- `bt` [`AbstractBetaFunction`](@ref) type struct
- `p_leaf` Leaf water potential `[MPa]`
- `p_soil` Soil water potential `[MPa]`
- `swc` Soil water content
"""
function β_factor(
            bt::Union{BetaGLinearPleaf{FT}, BetaVLinearPleaf{FT}},
            p_leaf::FT,
            p_soil::FT,
            swc::FT
) where {FT<:AbstractFloat}
    @unpack p_max, p_min = bt;

    if p_leaf >= p_max
        return FT(1)
    elseif p_leaf <= p_min
        return FT(0)
    else
        return (p_leaf - p_min) / (p_max - p_min)
    end
end




function β_factor(
            bt::Union{BetaGLinearPsoil{FT}, BetaVLinearPsoil{FT}},
            p_leaf::FT,
            p_soil::FT,
            swc::FT
) where {FT<:AbstractFloat}
    @unpack p_max, p_min = bt;

    if p_soil >= p_max
        return FT(1)
    elseif p_soil <= p_min
        return FT(0)
    else
        return (p_soil - p_min) / (p_max - p_min)
    end
end




function β_factor(
            bt::Union{BetaGLinearSWC{FT}, BetaVLinearSWC{FT}},
            p_leaf::FT,
            p_soil::FT,
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

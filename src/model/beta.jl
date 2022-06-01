###############################################################################
#
# Beta functions to make soil moisture correction
#
###############################################################################
"""
    β_factor(hs::LeafHydraulics{FT},
             svc::AbstractSoilVC{FT},
             bt::AbstractBetaFunction{FT},
             p_leaf::FT,
             p_soil::FT,
             swc::FT
    ) where {FT<:AbstractFloat}

Calculate the β correction factor, given
- `hs` `LeafHydraulics` structure
- `svc` Soil vulnerability curve
- `bt` [`AbstractBetaFunction`](@ref) type struct
- `p_leaf` Leaf xylem pressure corrected to 25 degree `[MPa]`
- `p_soil` Soil water potential corrected to 25 degree `[MPa]`
- `swc` Soil water content
"""
function β_factor(
            hs::LeafHydraulics{FT},
            svc::AbstractSoilVC{FT},
            bt::Union{BetaGLinearKleaf{FT}, BetaVLinearKleaf{FT}},
            p_leaf::FT,
            p_soil::FT,
            swc::FT
) where {FT<:AbstractFloat}
    @unpack VC = hs;

    return relative_hydraulic_conductance(VC, p_leaf)
end




function β_factor(
            hs::LeafHydraulics{FT},
            svc::AbstractSoilVC{FT},
            bt::Union{BetaGLinearKsoil{FT}, BetaVLinearKsoil{FT}},
            p_leaf::FT,
            p_soil::FT,
            swc::FT
) where {FT<:AbstractFloat}
    return relative_hydraulic_conductance(svc, true, p_soil)
end




function β_factor(
            hs::LeafHydraulics{FT},
            svc::AbstractSoilVC{FT},
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
            hs::LeafHydraulics{FT},
            svc::AbstractSoilVC{FT},
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
            hs::LeafHydraulics{FT},
            svc::AbstractSoilVC{FT},
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

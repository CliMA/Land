###############################################################################
#
# Temperature effects on surface tension and viscosity
#
###############################################################################
"""
    temperature_effects!(
                hs::LeafHydraulics{FT}
    ) where {FT<:AbstractFloat}
    temperature_effects!(
                hs::AbstractHydraulicOrgan{FT}
    ) where {FT<:AbstractFloat}
    temperature_effects!(
                hs::AbstractHydraulicOrgan{FT},
                T::FT
    ) where {FT<:AbstractFloat}
    temperature_effects!(
                tree::TreeSimple{FT}
    ) where {FT<:AbstractFloat}

Update temperature effetcs on hydralic system, given
- `hs` [`AbstractHydraulicOrgan`](@ref) type struct
- `T` Given temperature
- `tree` [`TreeSimple`](@ref) type struct
"""
function temperature_effects!(hs::LeafHydraulics{FT}) where {FT<:AbstractFloat}
    if hs.T_sap != hs.T_old
        hs.f_st  = relative_surface_tension(hs.T_sap);
        hs.f_vis = relative_viscosity(hs.T_sap);
        hs.p_crt = xylem_p_crit(hs.vc);
        hs.T_old = hs.T_sap;
    end

    return nothing
end




function temperature_effects!(hs::Union{RootHydraulics{FT}, StemHydraulics{FT}}) where {FT<:AbstractFloat}
    if hs.T_sap != hs.T_old
        hs.f_st  = relative_surface_tension(hs.T_sap);
        hs.f_vis = relative_viscosity(hs.T_sap);
        hs.T_old = hs.T_sap;
    end

    return nothing
end




function temperature_effects!(hs::Union{LeafHydraulics{FT}, RootHydraulics{FT}, StemHydraulics{FT}}, T::FT) where {FT<:AbstractFloat}
    if T != hs.T_sap
        hs.T_sap = T;
        temperature_effects!(hs);
    end

    return nothing
end




function temperature_effects!(tree::GrassLikeOrganism{FT}) where {FT<:AbstractFloat}
    for root in tree.roots
        temperature_effects!(root);
    end
    for leaf in tree.leaves
        temperature_effects!(leaf);
    end

    return nothing
end




function temperature_effects!(tree::PalmLikeOrganism{FT}) where {FT<:AbstractFloat}
    for root in tree.roots
        temperature_effects!(root);
    end
    temperature_effects!(tree.trunk);
    for leaf in tree.leaves
        temperature_effects!(leaf);
    end

    return nothing
end




function temperature_effects!(tree::TreeLikeOrganism{FT}) where {FT<:AbstractFloat}
    for root in tree.roots
        temperature_effects!(root);
    end
    temperature_effects!(tree.trunk);
    for stem in tree.branch
        temperature_effects!(stem);
    end
    for leaf in tree.leaves
        temperature_effects!(leaf);
    end

    return nothing
end




function temperature_effects!(tree::TreeSimple{FT}) where {FT<:AbstractFloat}
    temperature_effects!(tree.root);
    temperature_effects!(tree.stem);
    temperature_effects!(tree.leaf);

    return nothing
end

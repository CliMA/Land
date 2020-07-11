###############################################################################
#
# Temperature effects on surface tension and viscosity
#
###############################################################################
"""
    vc_temperature_effects!(hs::AbstractHydraulicSystem{FT}, T::FT) where {FT<:AbstractFloat}
    vc_temperature_effects!(hs::AbstractHydraulicSystem{FT}) where {FT<:AbstractFloat}
    vc_temperature_effects!(tree::TreeSimple{FT}) where {FT<:AbstractFloat}

Update temperature effetcs on hydralic system, given
- `hs` [`AbstractHydraulicSystem`](@ref) type struct
- `T` Given temperature
- `tree` [`TreeSimple`](@ref) type struct
"""
function vc_temperature_effects!(
            hs::LeafHydraulics{FT}
) where {FT<:AbstractFloat}
    if hs.T_sap != hs.T_old
        hs.f_st  = relative_surface_tension(hs.T_sap);
        hs.f_vis = relative_viscosity(hs.T_sap);
        hs.p_crt = xylem_p_crit(hs.vc);
        hs.T_old = hs.T_sap;
    end

    return nothing
end




function vc_temperature_effects!(
            hs::AbstractHydraulicSystem{FT}
) where {FT<:AbstractFloat}
    if hs.T_sap != hs.T_old
        hs.f_st  = relative_surface_tension(hs.T_sap);
        hs.f_vis = relative_viscosity(hs.T_sap);
        hs.T_old = hs.T_sap;
    end

    return nothing
end




function vc_temperature_effects!(
            hs::AbstractHydraulicSystem{FT},
            T::FT
) where {FT<:AbstractFloat}
    if T != hs.T_sap
        hs.T_sap = T;
        vc_temperature_effects!(hs);
    end

    return nothing
end




function vc_temperature_effects!(
            tree::TreeSimple{FT}
) where {FT<:AbstractFloat}
    vc_temperature_effects!(tree.root);
    vc_temperature_effects!(tree.stem);
    vc_temperature_effects!(tree.leaf);

    return nothing
end

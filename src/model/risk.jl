###############################################################################
#
# Calculate xylem end risk factor (k_ratio)
#
###############################################################################
"""
    xylem_risk(hs::LeafHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}

Evaluate the hydraulic risk at the end of leaf xylem, given
- `hs` `LeafHydraulics` type struct
- `flow` Flow rate (per leaf area)
- `T` Liquid temperature
"""
function xylem_risk(
            hs::LeafHydraulics{FT},
            flow::FT,
            T::FT
) where {FT<:AbstractFloat}
    @unpack VC = hs;

    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_25 = xylem_end_pressure(hs, flow, T);

    return relative_hydraulic_conductance(VC, _p_25);
end

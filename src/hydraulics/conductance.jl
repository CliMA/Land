###############################################################################
#
# Calculate xylem end risk factor (k_ratio)
#
###############################################################################
"""
    xylem_risk(hs::LeafHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}

Evaluate the hydraulic risk at the end of leaf xylem, given
- `hs` [`LeafHydraulics`](@ref) type struct
- `flow` Flow rate (per leaf area)
"""
function xylem_risk(
            hs::LeafHydraulics{FT},
            flow::FT
) where {FT<:AbstractFloat}
    @unpack f_st, f_vis, vc = hs;

    p_25 = end_pressure(hs, flow) / hs.f_st;
    k_25 = xylem_k_ratio(vc, p_25, f_vis);

    return k_25
end
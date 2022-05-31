#=
###############################################################################
#
# Whole plant hydraulics --- Element and whole-tree conductance ratio
#
###############################################################################
"""
    plant_conductances!(tree::TreeSimple{FT}) where {FT<:AbstractFloat}

Update plant total hydraulic information in each element and whole plant, given
- `tree` [`TreeSimple`](@ref) type struct
"""
function plant_conductances!(tree::MonoElementSPAC{FT}) where {FT<:AbstractFloat}
    # Calculate the resistances
    r_m_root = 1 / (tree.root).k_max;
    r_m_stem = 1 / (tree.stem).k_max;
    r_m_leaf = 1 / ((tree.leaf).k_sla * (tree.leaf).area);
    r_m_tree = r_m_root + r_m_stem + r_m_leaf

    # Calculate the relative resistances
    r_root = FT(0);
    r_stem = FT(0);
    r_leaf = FT(0);
    for i in 1:tree.root.N
        r_root += (tree.root).k_history[i];
    end
    for i in 1:tree.stem.N
        r_stem += (tree.stem).k_history[i];
    end
    for i in 1:tree.leaf.N
        r_leaf += (tree.leaf).k_history[i];
    end
    r_root /= tree.root.N;
    r_stem /= tree.stem.N;
    r_leaf /= tree.leaf.N;
    r_tree  = r_root * r_m_root + r_stem * r_m_stem + r_leaf * r_m_leaf;

    # Update the kr info
    tree.krs[1] = 1 / r_root;
    tree.krs[2] = 1 / r_stem;
    tree.krs[3] = 1 / r_leaf;
    tree.krs[4] = r_m_tree / r_tree;

    return nothing
end








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

    p_25 = xylem_end_pressure(hs, flow) / hs.f_st;
    T_25 = relative_hydraulic_conductance(vc, p_25) / f_vis;

    return T_25
end
=#

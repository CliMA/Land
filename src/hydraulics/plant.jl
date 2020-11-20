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
function plant_conductances!(
            tree::TreeSimple{FT}
) where {FT<:AbstractFloat}
    # Calculate the resistances
    r_m_root = 1 / (tree.root).k_max;
    r_m_stem = 1 / (tree.stem).k_max;
    r_m_leaf = 1 / ((tree.leaf).k_sla * (tree.leaf).area);
    r_m_tree = r_m_root + r_m_stem + r_m_leaf

    # Calculate the relative resistances
    r_root = FT(0);
    r_stem = FT(0);
    r_leaf = FT(0);
    for i in 1:10
        r_root += (tree.root).k_history[i];
        r_stem += (tree.stem).k_history[i];
        r_leaf += (tree.leaf).k_history[i];
    end
    r_root /= 10;
    r_stem /= 10;
    r_leaf /= 10;
    r_tree  = r_root * r_m_root + r_stem * r_m_stem + r_leaf * r_m_leaf;

    # Update the kr info
    tree.krs[1] = 1 / r_root;
    tree.krs[2] = 1 / r_stem;
    tree.krs[3] = 1 / r_leaf;
    tree.krs[4] = r_m_tree / r_tree;

    return nothing
end

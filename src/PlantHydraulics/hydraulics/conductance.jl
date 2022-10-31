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
function plant_conductances!(tree::TreeSimple{FT}) where {FT<:AbstractFloat}
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
# Root pressure and conductance
#
###############################################################################
"""
    root_pk(root::RootHydraulics{FT},
            flow::FT
    ) where {FT<:AbstractFloat}
    root_pk(root::RootHydraulics{FT},
            q_in::FT,
            flow::Array{FT,1}
    ) where {FT<:AbstractFloat}

Return root xylem end pressure and root hydraulic conductance (reverse of
    summed resistance), given
- `root` [`RootHydraulics`](@ref) struct
- `flow` Given flow rate(s) in the root layer, array for non-steady state with
    capacitance enabled
- `q_in` Flow rate into the root
"""
function root_pk(
            root::RootHydraulics{FT},
            flow::FT
) where {FT<:AbstractFloat}
    @unpack f_st, f_vis, k_element, k_history, k_rhiz, p_gravity, p_history,
            p_osm, p_ups, sh, T_sap, vc = root;

    # make sure that p_ups is not p_25 and then convert
    p_end::FT = p_ups;
    p_25 ::FT = p_end / f_st;
    r_all::FT = FT(0);

    # compute pressure drop along rhizosphere, using p_25 for Θ
    _dr = 1 / k_rhiz * f_vis / 10;
    _dp = flow * _dr;
    for i in 1:10
        # No idea why soil_k_ratio_p25 results in unnecessary allocations
        # _f   = soil_k_ratio_p25(sh, p_25);
        _rwc   = soil_rwc(sh, p_25);
        _f     = soil_k_ratio_rwc(sh, _rwc);
        p_25  -= _dp / _f;
        r_all += _dr / _f;
    end
    p_end = p_25 * f_st + p_osm * T_sap / T₂₅(FT);

    # compute k from temperature and history, then update pressure
    for (_k, _kh, _pg, _ph) in zip(k_element, k_history, p_gravity, p_history)
        p_25 = p_end / f_st;
        if p_25 < _ph
            k = xylem_k_ratio(vc, p_25, f_vis) * _k;
        else
            k = _kh * _k;
        end
        p_end -= flow / k + _pg;
        r_all += 1 / k;
    end

    k_all = 1 / r_all;

    return p_end, k_all
end




function root_pk(
            root::RootHydraulics{FT},
            q_in::FT,
            flow::Array{FT,1}
) where {FT<:AbstractFloat}
    @unpack f_st, f_vis, k_element, k_history, k_rhiz, p_gravity, p_history,
            p_osm, p_ups, sh, T_sap, vc = root;

    # make sure that p_ups is not p_25 and then convert
    p_end::FT = p_ups;
    p_25 ::FT = p_end / f_st;
    r_all::FT = FT(0);

    # compute pressure drop along rhizosphere, using p_25 for Θ
    _dr = 1 / k_rhiz * f_vis / 10;
    _dp = q_in * _dr;
    for i in 1:10
        # No idea why soil_k_ratio_p25 results in unnecessary allocations
        # _f   = soil_k_ratio_p25(sh, p_25);
        _rwc   = soil_rwc(sh, p_25);
        _f     = soil_k_ratio_rwc(sh, _rwc);
        p_25  -= _dp / _f;
        r_all += _dr / _f;
    end
    p_end = p_25 * f_st + p_osm * T_sap / T₂₅(FT);

    # compute k from temperature and history, then update pressure
    for (_k, _kh, _pg, _ph, _fl) in zip(k_element, k_history, p_gravity,
                                        p_history, flow)
        p_25 = p_end / f_st;
        if p_25 < _ph
            k = xylem_k_ratio(vc, p_25, f_vis) * _k;
        else
            k = _kh * _k;
        end
        p_end -= _fl / k + _pg;
        r_all += 1 / k;
    end

    k_all = 1 / r_all;

    return p_end, k_all
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

    p_25 = end_pressure(hs, flow) / hs.f_st;
    θ_25 = xylem_k_ratio(vc, p_25, f_vis);

    return θ_25
end

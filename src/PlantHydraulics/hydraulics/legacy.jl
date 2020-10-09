###############################################################################
#
# Update pressure profile along a hydraulic system
#
###############################################################################
"""
    hydraulic_p_profile!(leaf::LeafHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}
    hydraulic_p_profile!(root::RootHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}
    hydraulic_p_profile!(stem::StemHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}
    hydraulic_p_profile!(tree::TreeSimple{FT}, p_soil::FT, flow::FT) where {FT<:AbstractFloat}
    hydraulic_p_profile!(tree::TreeSimple{FT}, p_soil::FT, f_sl::FT, f_sh::FT, r_sl::FT) where {FT<:AbstractFloat}

Update the pressure profile, given
- `leaf` [`LeafHydraulics`](@ref) type struct
- `root` [`RootHydraulics`](@ref) type struct
- `stem` [`StemHydraulics`](@ref) type struct
- `tree` [`TreeSimple`](@ref) type struct
- `p_soil` Soil water potential
- `flow` Flow rate (per leaf area for [`LeafHydraulics`](@ref))
- `f_sl` Flow rate to sunlit leaves
- `f_sh` Flow rate to shaded leaves
- `r_sl` Fraction of sunlit leaves
"""
function hydraulic_p_profile!(
            leaf::LeafHydraulics{FT},
            flow::FT
) where {FT<:AbstractFloat}
    @unpack k_element, k_history, p_history, p_ups, f_st, f_vis, vc = leaf;

    p_end::FT = p_ups;

    # compute k from temperature and history, then update pressure
    for i in eachindex(k_element)
        # update history first
        p_25 = p_end / f_st;
        if p_25 < p_history[i]
            _kr = xylem_k_ratio(vc, p_25, f_vis);
            leaf.p_history[i] = p_25;
            leaf.k_history[i] = _kr;
            k = _kr * k_element[i];
        end

        # then use historical minimal k
        k = k_history[i] * k_element[i];
        p_end -= flow / k;

        leaf.p_element[i] = p_end;
    end

    # update the leaf xylem end pressure
    leaf.p_dos  = p_end;
    leaf.p_leaf = p_end - flow / leaf.k_ox;

    return nothing
end




function hydraulic_p_profile!(
            root::RootHydraulics{FT},
            flow::FT
) where {FT<:AbstractFloat}
    @unpack f_st, f_vis, k_element, k_history, k_rhiz, p_gravity, p_history,
            p_ups, sh, vc = root;

    # make sure that p_ups is not p_25 and then convert
    p_end::FT = p_ups;
    p_25 ::FT = p_end / f_st;

    # compute pressure drop along rhizosphere, using p_25 for Θ
    _dp = flow / k_rhiz * f_vis / 10;
    for i in 1:10
        # No idea why soil_k_ratio_p25 results in unnecessary allocations
        # _f  = soil_k_ratio_p25(sh, p_25);
        _rwc  = soil_rwc(sh, p_25);
        _f    = soil_k_ratio_rwc(sh, _rwc);
        p_25 -= _dp / _f;
    end
    p_end = p_25 * f_st;

    # compute k from temperature and history, then update pressure
    for i in eachindex(k_element)
        # update history first
        p_25 = p_end / f_st;
        if p_25 < p_history[i]
            _kr = xylem_k_ratio(vc, p_25, f_vis);
            root.p_history[i] = p_25;
            root.k_history[i] = _kr;
            k = _kr * k_element[i];
        end

        # then use historical minimal k
        k = k_history[i] * k_element[i];
        p_end -= flow / k + p_gravity[i];

        root.p_element[i] = p_end;
    end

    # update the leaf xylem end pressure
    root.p_dos  = p_end;

    return nothing
end




function hydraulic_p_profile!(
            stem::StemHydraulics{FT},
            flow::FT
) where {FT<:AbstractFloat}
    @unpack f_st, f_vis, k_element, k_history, p_gravity, p_history, p_ups,
            vc = stem;

    p_end::FT = p_ups;

    # compute k from temperature and history, then update pressure
    for i in eachindex(k_element)
        # update history first
        p_25 = p_end / f_st;
        if p_25 < p_history[i]
            _kr = xylem_k_ratio(vc, p_25, f_vis);
            stem.p_history[i] = p_25;
            stem.k_history[i] = _kr;
            k = _kr * k_element[i];
        end

        # then use historical minimal k
        k = k_history[i] * k_element[i];
        p_end -= flow / k + p_gravity[i];

        stem.p_element[i] = p_end;
    end

    # update the leaf xylem end pressure
    stem.p_dos  = p_end;

    return nothing
end




function hydraulic_p_profile!(
            tree::TreeSimple{FT},
            p_soil::FT,
            flow::FT
) where {FT<:AbstractFloat}
    # update the profile in root
    (tree.root).p_ups = p_soil;
    hydraulic_p_profile!(tree.root, flow);

    # update the profile in stem
    (tree.stem).p_ups = (tree.root).p_dos;
    hydraulic_p_profile!(tree.stem, flow);

    # update the profile in leaf
    (tree.leaf).p_ups = (tree.stem).p_dos;
    hydraulic_p_profile!(tree.leaf, flow);

    return nothing
end




function hydraulic_p_profile!(
            tree::TreeSimple{FT},
            p_soil::FT,
            f_sl::FT,
            f_sh::FT,
            r_sl::FT
) where {FT<:AbstractFloat}
    flow = f_sl + f_sh;

    # update the profile in root
    (tree.root).p_ups = p_soil;
    hydraulic_p_profile!(tree.root, flow);

    # update the profile in stem
    (tree.stem).p_ups = (tree.root).p_dos;
    hydraulic_p_profile!(tree.stem, flow);

    # update the profile in sunlit leaf
    (tree.leaf).p_ups = (tree.stem).p_dos;
    hydraulic_p_profile!(tree.leaf, f_sl / (tree.leaf).area / r_sl );
    hydraulic_p_profile!(tree.leaf, f_sh / (tree.leaf).area / (1-r_sl) );

    return nothing
end








###############################################################################
#
# Whole plant hydraulics --- Element and whole-tree conductance ratio
#
###############################################################################
"""
    inititialize_legacy!(hs::AbstractHydraulicSystem{FT}) where {FT<:AbstractFloat}
    inititialize_legacy!(tree::TreeSimple{FT}) where {FT<:AbstractFloat}

Initialize the drought legacy effects in the xylem, given
- `hs` [`AbstractHydraulicSystem`] type struct
- `tree` [`TreeSimple`](@ref) type struct
"""
function inititialize_legacy!(
            hs::AbstractHydraulicSystem{FT}
) where {FT<:AbstractFloat}
    hs.k_history .= 1
    hs.p_history .= 0

    return nothing
end




function inititialize_legacy!(tree::TreeSimple{FT}) where {FT<:AbstractFloat}
    inititialize_legacy!(tree.root);
    inititialize_legacy!(tree.stem);
    inititialize_legacy!(tree.leaf);

    return nothing
end

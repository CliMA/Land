###############################################################################
#
# Calculate pressure profile along a hydraulic system
#
###############################################################################
"""
    xylem_p_from_flow(leaf::LeafHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}
    xylem_p_from_flow(root::RootHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}
    xylem_p_from_flow(stem::StemHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}
    xylem_p_from_flow(tree::GrassLikeHS{FT}, flow::FT) where {FT<:AbstractFloat}
    xylem_p_from_flow(tree::PalmLikeHS{FT}, flow::FT) where {FT<:AbstractFloat}
    xylem_p_from_flow(tree::TreeLikeHS{FT}, flow::FT) where {FT<:AbstractFloat}
    xylem_p_from_flow(tree::TreeSimple{FT}, flow::FT) where {FT<:AbstractFloat}
    xylem_p_from_flow(tree::TreeSimple{FT}, f_sl::FT, f_sh::FT, r_sl::FT) where {FT<:AbstractFloat}

Return the xylen end pressure(s) from flow rate(s), given
- `leaf` [`LeafHydraulics`](@ref) type struct
- `root` [`RootHydraulics`](@ref) type struct
- `stem` [`StemHydraulics`](@ref) type struct
- `tree` [`AbstractPlantHS`](@ref) type struct
- `flow` Flow rate (per leaf area for [`LeafHydraulics`](@ref))
- `f_sl` Flow rate to sunlit leaves
- `f_sh` Flow rate to shaded leaves
- `r_sl` Fraction of sunlit leaves

Note, gravity is accounted for in root and stem; rhizosphere conductance is
    accounted for in root; extra-xylary conductance is not accounted for in
    leaf here because it calculates xylem end pressure.
"""
function xylem_p_from_flow(
            leaf::LeafHydraulics{FT},
            flow::FT
) where {FT<:AbstractFloat}
    @unpack f_st, f_vis, k_element, k_history, p_history, p_ups, vc = leaf;

    p_end::FT = p_ups;

    # compute k from temperature and history, then update pressure
    for (_k, _kh, _ph) in zip(k_element, k_history, p_history)
        p_25 = p_end / f_st;
        if p_25 < _ph
            k = xylem_k_ratio(vc, p_25, f_vis) * _k;
        else
            k = _kh * _k;
        end
        p_end -= flow / k;
    end

    return p_end
end




function xylem_p_from_flow(
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
    for (_k, _kh, _pg, _ph) in zip(k_element, k_history, p_gravity, p_history)
        p_25 = p_end / f_st;
        if p_25 < _ph
            k = xylem_k_ratio(vc, p_25, f_vis) * _k;
        else
            k = _kh * _k;
        end
        p_end -= flow / k + _pg;
    end

    return p_end
end




function xylem_p_from_flow(
            stem::StemHydraulics{FT},
            flow::FT
) where {FT<:AbstractFloat}
    @unpack f_st, f_vis, k_element, k_history, p_gravity, p_history, p_ups,
            vc = stem;

    p_end::FT = p_ups;

    # compute k from temperature and history, then update pressure
    for (_k, _kh, _pg, _ph) in zip(k_element, k_history, p_gravity, p_history)
        p_25 = p_end / f_st;
        if p_25 < _ph
            k = xylem_k_ratio(vc, p_25, f_vis) * _k;
        else
            k = _kh * _k;
        end
        p_end -= flow / k + _pg;
    end

    return p_end
end




function xylem_p_from_flow(
            tree::GrassLikeHS{FT},
            flow::FT
) where {FT<:AbstractFloat}
    # Total leaf area
    #tla = sum([tree.leaves[i].area for i in 1:tree.n_canopy])
    tla = FT(0)
    for leaf in tree.leaves
        tla += leaf.area;
    end

    if tree.n_canopy == 1
        # calculate the p_dos for roots
        roots_flow!(tree.roots, tree.container_k, tree.container_p, tree.container_q, flow);
        tree.leaves[1].p_ups = mean(tree.container_p);

        # calculate the p_dos for leaves
        p_dos = xylem_p_from_flow(tree.leaves[1], flow/tla);

        return p_dos
    else
        println("No function applicable for multi-layer canopy!")
        return ErrorException("Error!");
    end
end




function xylem_p_from_flow(
            tree::PalmLikeHS{FT},
            flow::FT
) where {FT<:AbstractFloat}
    # Total leaf area
    #tla = sum([tree.leaves[i].area for i in 1:tree.n_canopy])
    tla = FT(0)
    for leaf in tree.leaves
        tla += leaf.area;
    end

    if tree.n_canopy == 1
        # calculate the p_dos for roots
        roots_flow!(tree.roots, tree.container_k, tree.container_p, tree.container_q, flow);
        (tree.trunk).p_ups = mean(tree.container_p);

        # calculate the p_dos for trunk
        p_dos = xylem_p_from_flow(tree.trunk, flow);
        tree.leaves[1].p_ups = p_dos;

        # calculate the p_dos for leaves
        p_dos = xylem_p_from_flow(tree.leaves[1], flow/tla);

        return p_dos
    else
        println("No function applicable for multi-layer canopy!")
        return ErrorException("Error!");
    end
end




function xylem_p_from_flow(
            tree::TreeLikeHS{FT},
            flow::FT
) where {FT<:AbstractFloat}
    # Total leaf area
    #tla = sum([tree.leaves[i].area for i in 1:tree.n_canopy])
    tla = FT(0)
    for leaf in tree.leaves
        tla += leaf.area;
    end

    if tree.n_canopy == 1
        # calculate the p_dos for roots
        roots_flow!(tree.roots, tree.container_k, tree.container_p, tree.container_q, flow);
        (tree.trunk).p_ups = mean(tree.container_p);

        # calculate the p_dos for trunk
        p_dos = xylem_p_from_flow(tree.trunk, flow);
        tree.branch[1].p_ups = p_dos;

        # calculate the p_dos for branch
        p_dos = xylem_p_from_flow(tree.branch[1], flow);
        tree.leaves[1].p_ups = p_dos;

        # calculate the p_dos for leaves
        p_dos = xylem_p_from_flow(tree.leaves[1], flow/tla);

        return p_dos
    else
        println("No function applicable for multi-layer canopy!")
        return ErrorException("Error!");
    end
end




function xylem_p_from_flow(
            tree::TreeSimple{FT},
            flow::FT
) where {FT<:AbstractFloat}
    # calculate the p_dos for roots
    p_dos = xylem_p_from_flow(tree.root, flow);
    (tree.stem).p_ups = p_dos;

    # calculate the p_dos for stem
    p_dos = xylem_p_from_flow(tree.stem, flow);
    (tree.leaf).p_ups = p_dos;

    # calculate the p_dos for leaves
    p_dos = xylem_p_from_flow(tree.leaf, flow / (tree.leaf).area);

    return p_dos
end




function xylem_p_from_flow(
            tree::TreeSimple{FT},
            f_sl::FT,
            f_sh::FT,
            r_sl::FT
) where {FT<:AbstractFloat}
    flow = f_sl + f_sh;

    # calculate the p_dos for roots
    p_dos = xylem_p_from_flow(tree.root, flow);
    (tree.stem).p_ups = p_dos;

    # calculate the p_dos for stem
    p_dos = xylem_p_from_flow(tree.stem, flow);
    (tree.leaf).p_ups = p_dos;

    # calculate the p_dos for leaves
    p_sl = xylem_p_from_flow(tree.leaf, f_sl / (tree.leaf).area / r_sl);
    p_sh = xylem_p_from_flow(tree.leaf, f_sh / (tree.leaf).area / (1-r_sl));

    return p_sl,p_sh
end

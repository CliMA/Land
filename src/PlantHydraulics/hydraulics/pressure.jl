###############################################################################
#
# Calculate end pressure of a hydraulic system
#
###############################################################################
"""
    end_pressure(
                leaf::LeafHydraulics{FT},
                flow::FT
    ) where {FT<:AbstractFloat}
    end_pressure(
                root::RootHydraulics{FT},
                flow::FT
    ) where {FT<:AbstractFloat}
    end_pressure(
                stem::StemHydraulics{FT},
                flow::FT
    ) where {FT<:AbstractFloat}
    end_pressure(
                tree::TreeSimple{FT},
                flow::FT
    ) where {FT<:AbstractFloat}
    end_pressure(
                tree::TreeSimple{FT},
                f_sl::FT,
                f_sh::FT,
                r_sl::FT
    ) where {FT<:AbstractFloat}

Return the xylen end pressure(s) from flow rate(s), given
- `leaf` [`LeafHydraulics`](@ref) type struct
- `flow` Flow rate (per leaf area for [`LeafHydraulics`](@ref))
- `root` [`RootHydraulics`](@ref) type struct
- `stem` [`StemHydraulics`](@ref) type struct
- `tree` [`TreeSimple`](@ref) type struct
- `f_sl` Flow rate to sunlit leaves
- `f_sh` Flow rate to shaded leaves
- `r_sl` Fraction of sunlit leaves

Note, gravity is accounted for in root and stem; rhizosphere conductance is
    accounted for in root; extra-xylary conductance is not accounted for in
    leaf here because it calculates xylem end pressure.
"""
function end_pressure(
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
            k = _kh * _k / f_vis;
        end
        p_end -= flow / k;
    end

    return p_end
end




function end_pressure(
            root::RootHydraulics{FT},
            flow::FT
) where {FT<:AbstractFloat}
    @unpack f_st, f_vis, k_element, k_history, k_rhiz, p_gravity, p_history,
            p_osm, p_ups, sh, T_sap, vc = root;

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
    p_end = p_25 * f_st + p_osm * T_sap / T₂₅(FT);

    # compute k from temperature and history, then update pressure
    for (_k, _kh, _pg, _ph) in zip(k_element, k_history, p_gravity, p_history)
        p_25 = p_end / f_st;
        if p_25 < _ph
            k = xylem_k_ratio(vc, p_25, f_vis) * _k;
        else
            k = _kh * _k / f_vis;
        end
        p_end -= flow / k + _pg;
    end

    return p_end
end




function end_pressure(
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
            k = _kh * _k / f_vis;
        end
        p_end -= flow / k + _pg;
    end

    return p_end
end




function end_pressure(
            tree::TreeSimple{FT},
            flow::FT
) where {FT<:AbstractFloat}
    # calculate the p_dos for roots
    p_dos = end_pressure(tree.root, flow);
    (tree.stem).p_ups = p_dos;

    # calculate the p_dos for stem
    p_dos = end_pressure(tree.stem, flow);
    (tree.leaf).p_ups = p_dos;

    # calculate the p_dos for leaves
    p_dos = end_pressure(tree.leaf, flow / (tree.leaf).area);

    return p_dos
end




function end_pressure(
            tree::TreeSimple{FT},
            f_sl::FT,
            f_sh::FT,
            r_sl::FT
) where {FT<:AbstractFloat}
    flow = f_sl + f_sh;

    # calculate the p_dos for roots
    p_dos = end_pressure(tree.root, flow);
    (tree.stem).p_ups = p_dos;

    # calculate the p_dos for stem
    p_dos = end_pressure(tree.stem, flow);
    (tree.leaf).p_ups = p_dos;

    # calculate the p_dos for leaves
    p_sl = end_pressure(tree.leaf, f_sl / (tree.leaf).area / r_sl);
    p_sh = end_pressure(tree.leaf, f_sh / (tree.leaf).area / (1-r_sl));

    return p_sl,p_sh
end








###############################################################################
#
# Update pressure profile along a hydraulic system
#
###############################################################################
"""
    pressure_profile!(
                leaf::LeafHydraulics{FT},
                flow::FT;
                update::Bool = true
    ) where {FT<:AbstractFloat}
    pressure_profile!(
                root::RootHydraulics{FT},
                flow::FT;
                update::Bool = true
    ) where {FT<:AbstractFloat}
    pressure_profile!(
                root::RootHydraulics{FT},
                q_in::FT,
                flow::Array{FT,1};
                update::Bool = true
    ) where {FT<:AbstractFloat}
    pressure_profile!(
                stem::StemHydraulics{FT},
                flow::FT;
                update::Bool = true
    ) where {FT<:AbstractFloat}
    pressure_profile!(
                stem::StemHydraulics{FT},
                flow::Array{FT,1};
                update::Bool = true
    ) where {FT<:AbstractFloat}
    pressure_profile!(
                tree::TreeSimple{FT},
                p_soil::FT,
                flow::FT;
                update::Bool = true
    ) where {FT<:AbstractFloat}
    pressure_profile!(
                tree::TreeSimple{FT},
                p_soil::FT,
                f_sl::FT,
                f_sh::FT,
                r_sl::FT;
                update::Bool = true
    ) where {FT<:AbstractFloat}

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
- `update` Optional. If true, update drought legacy. Default is `true` for
    [`TreeSimple`](@ref) but `false` for others using capacitance
"""
function pressure_profile!(
            leaf::LeafHydraulics{FT},
            flow::FT;
            update::Bool = true
) where {FT<:AbstractFloat}
    # update leaf flow
    leaf.flow = flow;

    @unpack k_element, k_history, p_history, p_ups, f_st, f_vis, vc = leaf;

    p_end::FT = p_ups;

    # compute k from temperature and history, then update pressure
    for i in eachindex(k_element)
        # update history first
        p_25 = p_end / f_st;
        if p_25 < p_history[i]
            _kr = xylem_k_ratio(vc, p_25, f_vis);
            if update
                leaf.p_history[i] = p_25;
                leaf.k_history[i] = _kr * f_vis;
            end
            k = _kr * k_element[i];
        else
            k = k_history[i] / f_vis * k_element[i];
        end

        # then use historical minimal k
        p_end -= flow / k;

        leaf.p_element[i] = p_end;
    end

    # update the leaf xylem end pressure
    leaf.p_dos  = p_end;
    leaf.p_leaf = p_end - flow / leaf.k_ox;

    return nothing
end




function pressure_profile!(
            root::RootHydraulics{FT},
            flow::FT;
            update::Bool = true
) where {FT<:AbstractFloat}
    # update root flow rate
    root.flow = flow;

    @unpack f_st, f_vis, k_element, k_history, k_rhiz, p_gravity, p_history,
            p_osm, p_ups, sh, T_sap, vc = root;

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
    p_end = p_25 * f_st + p_osm * T_sap / T₂₅(FT);
    root.p_rhiz = p_end;

    # compute k from temperature and history, then update pressure
    for i in eachindex(k_element)
        # update history first
        p_25 = p_end / f_st;
        if p_25 < p_history[i]
            _kr = xylem_k_ratio(vc, p_25, f_vis);
            if update
                root.p_history[i] = p_25;
                root.k_history[i] = _kr * f_vis;
            end
            k = _kr * k_element[i];
        else
            k = k_history[i] / f_vis * k_element[i];
        end

        # then use historical minimal k
        p_end -= flow / k + p_gravity[i];

        root.p_element[i] = p_end;
    end

    # update the leaf xylem end pressure
    root.p_dos  = p_end;

    return nothing
end




function pressure_profile!(
            root::RootHydraulics{FT},
            q_in::FT,
            flow::Array{FT,1};
            update::Bool = true
) where {FT<:AbstractFloat}
    @unpack f_st, f_vis, k_element, k_history, k_rhiz, p_gravity, p_history,
            p_osm, p_ups, sh, T_sap, vc = root;

    # make sure that p_ups is not p_25 and then convert
    p_end::FT = p_ups;
    p_25 ::FT = p_end / f_st;

    # compute pressure drop along rhizosphere, using p_25 for Θ
    _dp = q_in / k_rhiz * f_vis / 10;
    for i in 1:10
        # No idea why soil_k_ratio_p25 results in unnecessary allocations
        # _f  = soil_k_ratio_p25(sh, p_25);
        _rwc  = soil_rwc(sh, p_25);
        _f    = soil_k_ratio_rwc(sh, _rwc);
        p_25 -= _dp / _f;
    end
    p_end = p_25 * f_st + p_osm * T_sap / T₂₅(FT);
    root.p_rhiz = p_end;

    # compute k from temperature and history, then update pressure
    for i in eachindex(k_element)
        # update history first
        p_25 = p_end / f_st;
        if p_25 < p_history[i]
            _kr = xylem_k_ratio(vc, p_25, f_vis);
            if update
                root.p_history[i] = p_25;
                root.k_history[i] = _kr * f_vis;
            end
            k = _kr * k_element[i];
        else
            k = k_history[i] / f_vis * k_element[i];
        end

        # then use historical minimal k
        p_end -= flow[i] / k + p_gravity[i];

        root.p_element[i] = p_end;
    end

    # update the leaf xylem end pressure
    root.p_dos  = p_end;

    return nothing
end




function pressure_profile!(
            stem::StemHydraulics{FT},
            flow::FT;
            update::Bool = true
) where {FT<:AbstractFloat}
    # update stem flow
    stem.flow = flow;

    @unpack f_st, f_vis, k_element, k_history, p_gravity, p_history, p_ups,
            vc = stem;

    p_end::FT = p_ups;

    # compute k from temperature and history, then update pressure
    for i in eachindex(k_element)
        # update history first
        p_25 = p_end / f_st;
        if p_25 < p_history[i]
            _kr = xylem_k_ratio(vc, p_25, f_vis);
            if update
                stem.p_history[i] = p_25;
                stem.k_history[i] = _kr * f_vis;
            end
            k = _kr * k_element[i];
        else
            k = k_history[i] / f_vis * k_element[i];
        end

        # then use historical minimal k
        p_end -= flow / k + p_gravity[i];

        stem.p_element[i] = p_end;
    end

    # update the leaf xylem end pressure
    stem.p_dos  = p_end;

    return nothing
end




function pressure_profile!(
            stem::StemHydraulics{FT},
            flow::Array{FT,1};
            update::Bool = true
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
            if update
                stem.p_history[i] = p_25;
                stem.k_history[i] = _kr * f_vis;
            end
            k = _kr * k_element[i];
        else
            k = k_history[i] / f_vis * k_element[i];
        end

        # then use historical minimal k
        p_end -= flow[i] / k + p_gravity[i];

        stem.p_element[i] = p_end;
    end

    # update the leaf xylem end pressure
    stem.p_dos  = p_end;

    return nothing
end




function pressure_profile!(
            tree::TreeSimple{FT},
            p_soil::FT,
            flow::FT;
            update::Bool = true
) where {FT<:AbstractFloat}
    # update the profile in root
    (tree.root).p_ups = p_soil;
    pressure_profile!(tree.root, flow; update=update);

    # update the profile in stem
    (tree.stem).p_ups = (tree.root).p_dos;
    pressure_profile!(tree.stem, flow; update=update);

    # update the profile in leaf
    (tree.leaf).p_ups = (tree.stem).p_dos;
    pressure_profile!(tree.leaf, flow; update=update);

    return nothing
end




function pressure_profile!(
            tree::TreeSimple{FT},
            p_soil::FT,
            f_sl::FT,
            f_sh::FT,
            r_sl::FT;
            update::Bool = true
) where {FT<:AbstractFloat}
    flow = f_sl + f_sh;

    # update the profile in root
    (tree.root).p_ups = p_soil;
    pressure_profile!(tree.root, flow; update=update);

    # update the profile in stem
    (tree.stem).p_ups = (tree.root).p_dos;
    pressure_profile!(tree.stem, flow; update=update);

    # update the profile in sunlit leaf
    (tree.leaf).p_ups = (tree.stem).p_dos;
    pressure_profile!(tree.leaf, f_sl / (tree.leaf).area / r_sl;
                         update=update);
    pressure_profile!(tree.leaf, f_sh / (tree.leaf).area / (1-r_sl);
                         update=update);

    return nothing
end




function pressure_profile!(
            tree::GrassLikeOrganism{FT},
            mode::NonSteadyStateMode;
            update::Bool = false
) where {FT<:AbstractFloat}
    leaves = tree.leaves;
    roots  = tree.roots;

    # update the profile in roots
    p_mean::FT = 0;
    for root in roots
        pressure_profile!(root, root.q_in, root.q_element; update=update);
        p_mean += root.p_dos;
    end
    p_mean /= length(roots);

    # update the profile in leaf
    for leaf in leaves
        leaf.p_ups = p_mean;
        pressure_profile!(leaf, leaf.q_in; update=update);
    end

    return nothing
end




function pressure_profile!(
            tree::GrassLikeOrganism{FT},
            mode::SteadyStateMode;
            update::Bool = false
) where {FT<:AbstractFloat}
    leaves = tree.leaves;
    roots  = tree.roots;

    # update the profile in roots
    p_mean::FT = 0;
    for root in roots
        pressure_profile!(root, root.flow; update=update);
        p_mean += root.p_dos;
    end
    p_mean /= length(roots);

    # update the profile in leaf
    for leaf in leaves
        leaf.p_ups = p_mean;
        pressure_profile!(leaf, leaf.flow; update=update);
    end

    return nothing
end




function pressure_profile!(
            tree::PalmLikeOrganism{FT},
            mode::NonSteadyStateMode;
            update::Bool = false
) where {FT<:AbstractFloat}
    leaves = tree.leaves;
    roots  = tree.roots;
    trunk  = tree.trunk;

    # update the profile in roots
    p_mean::FT = 0;
    for root in roots
        pressure_profile!(root, root.q_in, root.q_element; update=update);
        p_mean += root.p_dos;
    end
    p_mean /= length(roots);

    # update the profile in trunk
    trunk.p_ups = p_mean;
    pressure_profile!(trunk, trunk.q_element; update=update);

    # update the profile in leaf
    for leaf in leaves
        leaf.p_ups = trunk.p_dos;
        pressure_profile!(leaf, leaf.q_in; update=update);
    end

    return nothing
end




function pressure_profile!(
            tree::PalmLikeOrganism{FT},
            mode::SteadyStateMode;
            update::Bool = false
) where {FT<:AbstractFloat}
    leaves = tree.leaves;
    roots  = tree.roots;
    trunk  = tree.trunk;

    # update the profile in roots
    p_mean::FT = 0;
    for root in roots
        pressure_profile!(root, root.flow; update=update);
        p_mean += root.p_dos;
    end
    p_mean /= length(roots);

    # update the profile in trunk
    trunk.p_ups = p_mean;
    pressure_profile!(trunk, trunk.flow; update=update);

    # update the profile in leaf
    for leaf in leaves
        leaf.p_ups = trunk.p_dos;
        pressure_profile!(leaf, leaf.flow; update=update);
    end

    return nothing
end




function pressure_profile!(
            tree::TreeLikeOrganism{FT},
            mode::NonSteadyStateMode;
            update::Bool = false
) where {FT<:AbstractFloat}
    @unpack branch, leaves = tree;
    roots = tree.roots;
    trunk = tree.trunk;

    # update the profile in roots
    p_mean::FT = 0;
    for root in roots
        pressure_profile!(root, root.q_in, root.q_element; update=update);
        p_mean += root.p_dos;
    end
    p_mean /= length(roots);

    # update the profile in trunk
    trunk.p_ups = p_mean;
    pressure_profile!(trunk, trunk.q_element; update=update);

    # update the profile in leaf
    for i in eachindex(leaves)
        stem = branch[i];
        leaf = leaves[i];
        stem.p_ups = trunk.p_dos;
        pressure_profile!(stem, stem.q_element; update=update);
        leaf.p_ups = stem.p_dos;
        pressure_profile!(leaf, leaf.q_in; update=update);
    end

    return nothing
end




function pressure_profile!(
            tree::TreeLikeOrganism{FT},
            mode::SteadyStateMode;
            update::Bool = false
) where {FT<:AbstractFloat}
    @unpack branch, leaves = tree;
    roots = tree.roots;
    trunk = tree.trunk;

    # update the profile in roots
    p_mean::FT = 0;
    for root in roots
        pressure_profile!(root, root.flow; update=update);
        p_mean += root.p_dos;
    end
    p_mean /= length(roots);

    # update the profile in trunk
    trunk.p_ups = p_mean;
    pressure_profile!(trunk, trunk.flow; update=update);

    # update the profile in leaf
    for i in eachindex(leaves)
        stem = branch[i];
        leaf = leaves[i];
        stem.p_ups = trunk.p_dos;
        pressure_profile!(stem, stem.flow; update=update);
        leaf.p_ups = stem.p_dos;
        pressure_profile!(leaf, leaf.flow; update=update);
    end

    return nothing
end

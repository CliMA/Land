###############################################################################
#
# Update pressure profile along a hydraulic system
#
###############################################################################
"""
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
        _f    = relative_hydraulic_conductance(sh, true, p_25);
        p_25 -= _dp / _f;
    end
    p_end = p_25 * f_st + p_osm * T_sap / T_25(FT);
    root.p_rhiz = p_end;

    # compute k from temperature and history, then update pressure
    for i in eachindex(k_element)
        # update history first
        p_25 = p_end / f_st;
        if p_25 < p_history[i]
            _kr = relative_hydraulic_conductance(vc, p_25) / f_vis;
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
            _kr = relative_hydraulic_conductance(vc, p_25) / f_vis;
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
            tree::MonoGrassSPAC{FT},
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
            tree::MonoPalmSPAC{FT},
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
            tree::MonoTreeSPAC{FT},
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

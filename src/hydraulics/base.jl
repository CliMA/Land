###############################################################################
#
# Calculate pressure profile along a hydraulic system
#
###############################################################################
"""
    xylem_p_from_flow(leaf::LeafHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}
    xylem_p_from_flow(root::RootHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}
    xylem_p_from_flow(stem::StemHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}

Return the xylen end pressure from flow rate, given
- `leaf` [`LeafHydraulics`](@ref) type struct
- `root` [`RootHydraulics`](@ref) type struct
- `stem` [`StemHydraulics`](@ref) type struct
- `flow` Flow rate (per leaf area)

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
        _f    = soil_k_ratio_p25(sh, p_25);
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








###############################################################################
#
# Update pressure profile along a hydraulic system
#
###############################################################################
"""
    hydraulic_p_profile!(leaf::LeafHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}
    hydraulic_p_profile!(root::RootHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}
    hydraulic_p_profile!(stem::StemHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}

Update the pressure profile, given
- `leaf` [`LeafHydraulics`](@ref) type struct
- `root` [`RootHydraulics`](@ref) type struct
- `stem` [`StemHydraulics`](@ref) type struct
- `flow` Flow rate (per leaf area)
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
        _f    = soil_k_ratio_p25(sh, p_25);
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

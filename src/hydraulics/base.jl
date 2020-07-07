###############################################################################
#
# Calculate pressure profile along a hydraulic system
#
###############################################################################
"""
    xylem_p_from_flow(hs::LeafHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}
    xylem_p_from_flow(hs::RootHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}
    xylem_p_from_flow(hs::StemHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}

Return the xylen end pressure from flow rate, given
- `hs` [`LeafHydraulics`](@ref) or [`RootHydraulics`](@ref) or
    [`StemHydraulics`](@ref) type struct
- `flow` Flow rate (per leaf area)

Note, gravity is accounted for in root and stem; rhizosphere conductance is
    accounted for in root; extra-xylary conductance is not accounted for in
    leaf here because it calculates xylem end pressure.
"""
function xylem_p_from_flow(
            hs::LeafHydraulics{FT},
            flow::FT
) where {FT<:AbstractFloat}
    @unpack f_st, f_vis, k_element, k_history, p_history, p_ups, vc = hs;

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
            hs::RootHydraulics{FT},
            flow::FT
) where {FT<:AbstractFloat}
    @unpack f_st, f_vis, k_element, k_history, k_rhiz, p_gravity, p_history,
            p_ups, soil_α, soil_m, soil_n, vc = hs;

    # make sure that p_ups is not p_25 and then convert
    p_end::FT = p_ups;
    p_25 ::FT = p_end / f_st;

    # compute pressure drop along rhizosphere, using p_25 for Θ
    _dp = flow / k_rhiz * f_vis / 10;
    for i in 1:10
        if p_25<=0
            _Θ = (1 / (1 + (soil_α*(-p_25))^soil_n)) ^ soil_m;
            _f = max(FT(1e-20), sqrt(_Θ) * (1 - (1-_Θ^(1/soil_m)) ^ soil_m)^2);
        else
            _f = FT(1);
        end
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
            hs::StemHydraulics{FT},
            flow::FT
) where {FT<:AbstractFloat}
    @unpack f_st, f_vis, k_element, k_history, p_gravity, p_history, p_ups,
            vc = hs;

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
    hydraulic_p_profile!(hs::LeafHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}
    hydraulic_p_profile!(hs::RootHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}
    hydraulic_p_profile!(hs::StemHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}

Update the pressure profile, given
- `hs` [`LeafHydraulics`](@ref) or [`RootHydraulics`](@ref) or
    [`StemHydraulics`](@ref) type struct
- `flow` Flow rate (per leaf area)
"""
function hydraulic_p_profile!(
            hs::LeafHydraulics{FT},
            flow::FT
) where {FT<:AbstractFloat}
    @unpack k_element, k_history, p_history, p_ups, f_st, f_vis, vc = hs;

    p_end::FT = p_ups;

    # compute k from temperature and history, then update pressure
    for i in eachindex(k_element)
        # update history first
        p_25 = p_end / f_st;
        if p_25 < p_history[i]
            _kr = xylem_k_ratio(vc, p_25, f_vis);
            hs.p_history[i] = p_25;
            hs.k_history[i] = _kr;
            k = _kr * k_element[i];
        end

        # then use historical minimal k
        k = k_history[i] * k_element[i];
        p_end -= flow / k;

        hs.p_element[i] = p_end;
    end

    # update the leaf xylem end pressure
    hs.p_dos  = p_end;
    hs.p_leaf = p_end - flow / hs.k_ox;

    return nothing
end




function hydraulic_p_profile!(
            hs::RootHydraulics{FT},
            flow::FT
) where {FT<:AbstractFloat}
    @unpack f_st, f_vis, k_element, k_history, k_rhiz, p_gravity, p_history,
            p_ups, soil_α, soil_m, soil_n, vc = hs;

    # make sure that p_ups is not p_25 and then convert
    p_end::FT = p_ups;
    p_25 ::FT = p_end / f_st;

    # compute pressure drop along rhizosphere, using p_25 for Θ
    _dp = flow / k_rhiz * f_vis / 10;
    for i in 1:10
        if p_25<=0
            _Θ = (1 / (1 + (soil_α*(-p_25))^soil_n)) ^ soil_m;
            _f = max(FT(1e-20), sqrt(_Θ) * (1 - (1-_Θ^(1/soil_m)) ^ soil_m)^2);
        else
            _f = FT(1);
        end
        p_25 -= _dp / _f;
    end
    p_end = p_25 * f_st;

    # compute k from temperature and history, then update pressure
    for i in eachindex(k_element)
        # update history first
        p_25 = p_end / f_st;
        if p_25 < p_history[i]
            _kr = xylem_k_ratio(vc, p_25, f_vis);
            hs.p_history[i] = p_25;
            hs.k_history[i] = _kr;
            k = _kr * k_element[i];
        end

        # then use historical minimal k
        k = k_history[i] * k_element[i];
        p_end -= flow / k + p_gravity[i];

        hs.p_element[i] = p_end;
    end

    # update the leaf xylem end pressure
    hs.p_dos  = p_end;

    return nothing
end




function hydraulic_p_profile!(
            hs::StemHydraulics{FT},
            flow::FT
) where {FT<:AbstractFloat}
    @unpack f_st, f_vis, k_element, k_history, p_gravity, p_history, p_ups,
            vc = hs;

    p_end::FT = p_ups;

    # compute k from temperature and history, then update pressure
    for i in eachindex(k_element)
        # update history first
        p_25 = p_end / f_st;
        if p_25 < p_history[i]
            _kr = xylem_k_ratio(vc, p_25, f_vis);
            hs.p_history[i] = p_25;
            hs.k_history[i] = _kr;
            k = _kr * k_element[i];
        end

        # then use historical minimal k
        k = k_history[i] * k_element[i];
        p_end -= flow / k + p_gravity[i];

        hs.p_element[i] = p_end;
    end

    # update the leaf xylem end pressure
    hs.p_dos  = p_end;

    return nothing
end

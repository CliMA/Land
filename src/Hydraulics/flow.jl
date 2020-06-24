###############################################################################
#
# Calculate pressure profile along a hydraulic system
#
###############################################################################
"""
    xylem_p_from_flow(hs::LeafHydraulics{FT}, flow::FT)

Return the xylen end pressure from flow rate, given
- `hs` [`AbstractHydraulicSystem`](@ref) type struct
- `flow` Flow rate (per leaf area)

Note, gravity is accounted for in root and stem; rhizosphere conductance is
    accounted for in root; extra-xylary conductance is not accounted for in
    leaf here because it calculates xylem end pressure.
"""
function xylem_p_from_flow(
            hs::LeafHydraulics{FT},
            flow::FT
            ) where {FT<:AbstractFloat}
    @unpack b, c, k_element, k_history, p_history, p_ups, f_st, f_vis = hs;

    p_end::FT = p_ups;

    # compute k from temperature and history, then update pressure
    for (_k, _kh, _ph) in zip(k_element, k_history, p_history)
        p_25 = p_end / f_st;
        if p_25 < _ph
            k = weibull_k_ratio(b, c, p_25, f_vis) * _k;
        else
            k = _kh * _k;
        end
        p_end -= flow / k;
    end

    return p_end
end








###############################################################################
#
# Update pressure profile along a hydraulic system
#
###############################################################################
"""
    hydraulic_p_profile!(hs::LeafHydraulics{FT}, flow::FT)

Update the pressure profile, given
- `hs` [`AbstractHydraulicSystem`](@ref) type struct
- `flow` Flow rate (per leaf area)
"""
function hydraulic_p_profile!(
            hs::LeafHydraulics{FT},
            flow::FT
            ) where {FT<:AbstractFloat}
    @unpack b, c, k_element, k_history, p_history, p_ups, f_st, f_vis = hs;

    p_end::FT = p_ups;

    # compute k from temperature and history, then update pressure
    for i in eachindex(k_element)
        # update history first
        p_25 = p_end / f_st;
        if p_25 < p_history[i]
            _kr = weibull_k_ratio(b, c, p_25, f_vis);
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








###############################################################################
#
# Calculate xylem end risk factor (k_ratio)
#
###############################################################################
"""
    leaf_xylem_risk(hs::LeafHydraulics{FT}, flow::FT)

Evaluate the hydraulic risk at the end of leaf xylem, given
- `hs` [`AbstractHydraulicSystem`](@ref) type struct
- `flow` Flow rate (per leaf area)
"""
function leaf_xylem_risk(
            hs::LeafHydraulics{FT},
            flow::FT
            ) where {FT<:AbstractFloat}
    @unpack b, c, f_st, f_vis = hs;

    p_25 = xylem_p_from_flow(hs, flow) / hs.f_st;
    k_25 = weibull_k_ratio(b, c, p_25, f_vis);

    return k_25
end








###############################################################################
#
# Calculate leaf-level e_crit
#
###############################################################################
"""
    leaf_e_crit(hs::LeafHydraulics{FT}, ini::FT = FT(2e-9))

Calculate the critical flow rate (K ≈ 0), given
- `hs` [`AbstractHydraulicSystem`](@ref) type struct
- `ini` Initial guess of `e_crit`, use last `e_crit` when possible

Note, for the safety of no NaN, update e_crit when ΔP >= -0.01
"""
function leaf_e_crit(
            hs::LeafHydraulics{FT},
            ini::FT = FT(2e-9)
            ) where {FT<:AbstractFloat}
    _fl = max(ini-FT(1e-6), FT(1e-9));
    _fh = max(ini     , FT(2e-9));
    _sm = SecantMethod{FT}(_fl, _fh);
    _cs = CompactSolution();
    _st = SolutionTolerance{FT}(1e-10);

    # use 1/p to make the curve concave
    @inline f(x) = 1/xylem_p_from_flow(hs, x) - 1/hs.p_crt;
    _solut  = find_zero(f, _sm, _cs, _st);
    _ec::FT = _solut.root

    return _ec
end

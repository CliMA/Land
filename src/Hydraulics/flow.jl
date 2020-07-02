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
    @unpack b, c, f_st, f_vis, k_element, k_history, p_history, p_ups = hs;

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

function xylem_p_from_flow(
            hs::RootHydraulics{FT},
            flow::FT
            ) where {FT<:AbstractFloat}
    @unpack b, c, f_st, f_vis, k_element, k_history, k_rhiz, p_gravity,
            p_history, p_ups, soil_α, soil_m, soil_n = hs;

    # TODO make sure that p_ups is not p_25 and then convert
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
            k = weibull_k_ratio(b, c, p_25, f_vis) * _k;
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
    @unpack b, c, f_st, f_vis, k_element, k_history, p_gravity, p_history,
            p_ups = hs;

    p_end::FT = p_ups;

    # compute k from temperature and history, then update pressure
    for (_k, _kh, _pg, _ph) in zip(k_element, k_history, p_gravity, p_history)
        p_25 = p_end / f_st;
        if p_25 < _ph
            k = weibull_k_ratio(b, c, p_25, f_vis) * _k;
        else
            k = _kh * _k;
        end
        p_end -= flow / k + _pg;
    end

    return p_end
end








###############################################################################
#
# Calculate root flow rates and pressure from total flow rate
#
###############################################################################
"""
    root_q_from_pressure(root::RootHydraulics, pressure::FT, ini::FT)

Calculate the flow rate from a given tree base pressure, given
- `root` [`RootHydraulics`](@ref) type struct
- `pressure` Given tree base pressure in `[MPa]`
- `ini` Start point for RootSolvers
"""
function root_q_from_pressure(
            root::RootHydraulics,
            pressure::FT,
            ini::FT = FT(2e-6)
            ) where {FT<:AbstractFloat}
    _fl = max(ini-1e-4, 1e-6);
    _fh = max(ini     , 2e-6);
    _sm = SecantMethod{FT}(_fl, _fh);
    _cs = CompactSolution();
    _st = SolutionTolerance{FT}(1e-7);

    # use -20 to avoid the values goes to negative
    @inline f(x) = max(xylem_p_from_flow(root, x), -20) - pressure;
    _solut = find_zero(f, _sm, _cs, _st);
    _q::FT = _solut.root;

    return _q
end




"""
    q_diff(roots::Array{RootHydraulics{FT},1}, pressure::FT, target::FT)

Function to calculate the flow rate difference to find solution pressure, given
- `root` [`RootHydraulics`](@ref) type struct
- `pressure` Given tree base pressure in `[MPa]`
- `target` Target flow rate in `[mol s⁻¹]`
"""
function q_diff(
            roots::Array{RootHydraulics{FT},1},
            pressure::FT,
            target::FT
            ) where {FT<:AbstractFloat}
    # initial guess is the stored flow rate
    _fs = [root.flow for root in roots]
    _qs = root_q_from_pressure.(roots, pressure, _fs);

    return sum(_qs) - target
end




"""
    root_qs_p_from_q(roots::Array{RootHydraulics{FT},1}, flow::FT, ini::FT)

Solve for the flow rates in each root layer and root pressure, given
- `roots` Array of [`RootHydraulics`](@ref) struts
- `flow` Total flow rates
- `ini` Initial guess of plant base pressure in `[MPa]`
"""
function root_qs_p_from_q(
            roots::Array{RootHydraulics{FT},1},
            flow::FT,
            ini::FT = FT(-2e-5)
            ) where {FT<:AbstractFloat}
    _pl = min(ini     , -2e-5);
    _ph = min(ini+1e-4, -1e-5);
    _sm = SecantMethod{FT}(_pl, _ph);
    _cs = CompactSolution();
    _st = SolutionTolerance{FT}(1e-6);

    @inline f(x) = q_diff(roots, x, flow);
    _solut = find_zero(f, _sm, _cs, _st);
    _p::FT = _solut.root;

    _qs = FT[root_q_from_pressure(root, _p, root.flow) for root in roots];

    return _qs, _p
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

function hydraulic_p_profile!(
            hs::RootHydraulics{FT},
            flow::FT
            ) where {FT<:AbstractFloat}
    @unpack b, c, f_st, f_vis, k_element, k_history, k_rhiz, p_gravity,
            p_history, p_ups, soil_α, soil_m, soil_n = hs;

    # TODO make sure that p_ups is not p_25 and then convert
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
            _kr = weibull_k_ratio(b, c, p_25, f_vis);
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
    @unpack b, c, f_st, f_vis, k_element, k_history, p_gravity, p_history,
            p_ups = hs;

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
        p_end -= flow / k + p_gravity[i];

        hs.p_element[i] = p_end;
    end

    # update the leaf xylem end pressure
    hs.p_dos  = p_end;

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
    leaf_e_crit(hs::LeafHydraulics{FT})

Calculate the critical flow rate (K ≈ 0), given
- `hs` [`AbstractHydraulicSystem`](@ref) type struct

Note, for the safety of no NaN, update e_crit when ΔP >= -0.01
"""
function leaf_e_crit(
            hs::LeafHydraulics{FT}
            ) where {FT<:AbstractFloat}
    _fl = FT(0);
    _fh = FT(10);
    _sm = NewtonBisectionMethod{FT}(_fl, _fh);
    _rt = ResidualTolerance{FT}(1e-5);

    # use hs.p_crt-1 to avoid the pressure being too negative
    @inline f(x) = xylem_p_from_flow(hs, x) - hs.p_crt;
    _solut  = find_zero_ext(f, _sm, _rt);
    _ec::FT = _solut;

    if isnan(_ec)
        println("E_crit is NaN, please check the settings...")
        @show hs.p_ups;
        println("");
    end

    return _ec
end

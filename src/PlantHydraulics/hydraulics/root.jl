###############################################################################
#
# Calculate root flow rates and pressure from total flow rate
#
###############################################################################
"""
    root_q_from_pressure(root::RootHydraulics{FT}, pressure::FT, ini::FT) where {FT<:AbstractFloat}

Calculate the flow rate from a given tree base pressure, given
- `root` [`RootHydraulics`](@ref) type struct
- `pressure` Given tree base pressure in `[MPa]`
- `ini` Initial guess
"""
function root_q_from_pressure(
            root::RootHydraulics{FT},
            pressure::FT,
            ini::FT = FT(1)
) where {FT<:AbstractFloat}
    _fh    = (root.p_ups - pressure) * root.k_max / root.f_vis;
    _fl    = -_fh;
    _fx    = min(_fh, ini);
    _ms    = NewtonBisectionMethod{FT}(_fl, _fh, _fx);
    _st    = SolutionTolerance{FT}(1e-4, 50);
    @inline f(x) = xylem_p_from_flow(root, x) - pressure;
    _solut = find_zero(f, _ms, _st);

    return _solut
end




"""
    root_pk_from_flow(root::RootHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}

Return root xylem end pressure and root hydraulic conductance (reverse of
    summed resistance), given
- `root` [`RootHydraulics`](@ref) struct
- `flow` Given flow rates in the root layer
"""
function root_pk_from_flow(
            root::RootHydraulics{FT},
            flow::FT
) where {FT<:AbstractFloat}
    @unpack f_st, f_vis, k_element, k_history, k_rhiz, p_gravity, p_history,
            p_ups, sh, vc = root;

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
        r_all += 1 / k;
    end

    k_all = 1 / r_all;

    return p_end, k_all
end




"""
    rrecalculate_roots_flow!(roots::Array{RootHydraulics{FT},1}, ks::Array{FT,1}, ps::Array{FT,1}, qs::Array{FT,1}, flow::FT) where {FT<:AbstractFloat}

Recalculate the flow rates in the root from the pressure and conductance
    profiles in each root at non-steady state, given
- `roots` Array of [`RootHydraulics`](@ref) structs
- `ks` Container for conductance in each root layer
- `ps` Container for end xylem pressure in each layer
- `qs` Container for flow rate in each layer
- `flow` Total flow rate
"""
function recalculate_roots_flow!(
            roots::Array{RootHydraulics{FT},1},
            ks::Array{FT,1},
            ps::Array{FT,1},
            qs::Array{FT,1},
            flow::FT
) where {FT<:AbstractFloat}
    N = length(roots);

    # flush the values to ks, ps, and qs
    for i in 1:N
        _p,_k = root_pk_from_flow(roots[i], qs[i]);
        ps[i] = _p;
        ks[i] = _k;
    end

    # use ps and ks to compute the Δq to adjust
    pm = mean(ps);
    for i in 1:N
        qs[i] -= (pm - ps[i]) * ks[i]
    end

    # adjust the qs so that sum(qs) = flow
    q_diff = sum(qs) - flow;
    k_sum  = sum(ks);
    for i in 1:N
        qs[i] -= q_diff * ks[1] / k_sum;
    end

    return nothing
end




"""
    root_qs_p_from_q(roots::Array{RootHydraulics{FT},1}, ks::Array{FT,1}, ps::Array{FT,1}, qs::Array{FT,1}, flow::FT) where {FT<:AbstractFloat}

Solve for the flow rates in each root layer and root pressure at steady state,
    given
- `roots` Array of [`RootHydraulics`](@ref) structs
- `ks` Container for conductance in each root layer
- `ps` Container for end xylem pressure in each layer
- `qs` Container for flow rate in each layer
- `flow` Total flow rate
"""
function roots_flow!(
            roots::Array{RootHydraulics{FT},1},
            ks::Array{FT,1},
            ps::Array{FT,1},
            qs::Array{FT,1},
            flow::FT
) where {FT<:AbstractFloat}
    count = 0
    while true
        recalculate_roots_flow!(roots, ks, ps, qs, flow);
        count += 1;

        if maximum(ps) - minimum(ps) < 1e-4
            break
        end

        if count>20
            break
        end
    end

    return nothing
end




function roots_flow!(
            plant::GrassLikeHS{FT},
            flow::FT
) where {FT<:AbstractFloat}
    roots_flow!(plant.roots,
                plant.container_k,
                plant.container_p,
                plant.container_q,
                flow);

    return nothing
end




function roots_flow!(
            plant::PalmLikeHS{FT},
            flow::FT
) where {FT<:AbstractFloat}
    roots_flow!(plant.roots,
                plant.container_k,
                plant.container_p,
                plant.container_q,
                flow);

    return nothing
end




function roots_flow!(
            plant::TreeLikeHS{FT},
            flow::FT
) where {FT<:AbstractFloat}
    roots_flow!(plant.roots,
                plant.container_k,
                plant.container_p,
                plant.container_q,
                flow);

    return nothing
end

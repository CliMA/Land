###############################################################################
#
# Capacitance buffer flow
#
###############################################################################
"""
    buffer_rate(pv::AbstractCapacity{FT}) where {FT<:AbstractFloat}

Return the buffer rate, given
- `pv` [`AbstractCapacity`](@ref) type struct

Note that only symplastic water can be used as capacitance
"""
function buffer_rate(
            pv::PVCurveLinear{FT}
) where {FT<:AbstractFloat}
    return pv.k_refill
end




function buffer_rate(
            pv::PVCurveSegmented{FT}
) where {FT<:AbstractFloat}
    return pv.k_refill * (1 - pv.RWC_apo)
end








###############################################################################
#
# Calculate leaf-level critical_flow
#
###############################################################################
"""
    critical_flow(
                hs::LeafHydraulics{FT},
                ini::FT = FT(0.5)
    ) where {FT<:AbstractFloat}
    critical_flow(
                tree::TreeSimple{FT},
                ini::FT = FT(0.5)
    ) where {FT<:AbstractFloat}

Calculate the critical flow rate (K ≈ 0), given
- `hs` [`LeafHydraulics`](@ref) type struct
- `ini` Initial guess
- `tree` [`TreeSimple`](@ref) type struct

Note, for the safety of no NaN, update critical_flow when ΔP >= -0.01
"""
function critical_flow(
            hs::LeafHydraulics{FT},
            ini::FT = FT(0.5)
) where {FT<:AbstractFloat}
    # calculate maximal flow
    _fh = (hs.p_ups - hs.p_crt) * hs.k_sla / hs.f_vis;
    _fl = FT(0);
    _fx = min((_fh+_fl)/2, ini);
    _ms = NewtonBisectionMethod{FT}(x_min=_fl, x_max=_fh, x_ini=_fx);
    _rt = SolutionTolerance{FT}(eps(FT)*100, 50);
    @inline f(x) = end_pressure(hs, x) - hs.p_crt;
    _solut  = find_zero(f, _ms, _rt);

    if isnan(_solut)
        @warn twarn("E_crit is NaN, please check the settings...") hs.p_ups;
    end

    return _solut
end




function critical_flow(
            tree::TreeSimple{FT},
            ini::FT = FT(0.5)
) where {FT<:AbstractFloat}
    # calculate maximal flow for whole tree, remember to times leaf area
    _kr = (tree.root).k_max;
    _ks = (tree.stem).k_max;
    _kl = (tree.leaf).k_sla / (tree.leaf).f_vis * (tree.leaf).area;
    _kt = 1 / (1 / _kr + 1 / _ks + 1/ _kl);
    _fh = -1 * (tree.leaf).p_crt * _kt;
    _fl = FT(0);
    _fx = min((_fh+_fl)/2, ini);
    _ms = NewtonBisectionMethod{FT}(x_min=_fl, x_max=_fh, x_ini=_fx);
    _rt = SolutionTolerance{FT}(eps(FT)*100, 50);
    @inline f(x) = end_pressure(tree, x) - (tree.leaf).p_crt;
    _solut  = find_zero(f, _ms, _rt);

    return _solut
end








###############################################################################
#
# Calculate root flow rates and pressure from total flow rate
#
###############################################################################
"""
    roots_flow!(roots::Array{RootHydraulics{FT},1},
                ks::Array{FT,1},
                ps::Array{FT,1},
                qs::Array{FT,1},
                flow::FT,
                recalculate::Bool
    ) where {FT<:AbstractFloat}
    roots_flow!(roots::Array{RootHydraulics{FT},1},
                ks::Array{FT,1},
                ps::Array{FT,1},
                qs::Array{FT,1},
                flow::FT
    ) where {FT<:AbstractFloat}
    roots_flow!(plant::Union{GrassLikeOrganism{FT},
                             PalmLikeOrganism{FT},
                             TreeLikeOrganism{FT}},
                flow::FT
    ) where {FT<:AbstractFloat}

Recalculate the flow rates in the root from the pressure and conductance
    profiles in each root at non-steady state, given
- `roots` Array of [`RootHydraulics`](@ref) structs
- `ks` Container for conductance in each root layer
- `ps` Container for end xylem pressure in each layer
- `qs` Container for flow rate out of each layer
- `flow` Total flow rate out of the roots
- `recalculate` A paceholder indicator of recalculating root flow (useless)
- `plant` [`AbstractPlantOrganism`](@ref) type struct
"""
function roots_flow!(
            roots::Array{RootHydraulics{FT},1},
            ks::Array{FT,1},
            ps::Array{FT,1},
            qs::Array{FT,1},
            flow::FT,
            recalculate::Bool
) where {FT<:AbstractFloat}
    N = length(roots);

    # flush the values to ks, ps, and qs
    for i in 1:N
        _p,_k = root_pk(roots[i], qs[i]);
        ps[i] = _p;
        ks[i] = _k;
    end

    # use ps and ks to compute the Δq to adjust
    pm = mean(ps);
    for i in 1:N
        qs[i] -= (pm - ps[i]) * ks[i];
    end

    # adjust the qs so that sum(qs) = flow
    q_diff = sum(qs) - flow;
    k_sum  = sum(ks);
    for i in 1:N
        qs[i] -= q_diff * ks[1] / k_sum;
    end

    return nothing
end




function roots_flow!(
            roots::Array{RootHydraulics{FT},1},
            ks::Array{FT,1},
            ps::Array{FT,1},
            qs::Array{FT,1},
            flow::FT
) where {FT<:AbstractFloat}
    count = 0
    while true
        roots_flow!(roots, ks, ps, qs, flow, true);
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
            plant::Union{GrassLikeOrganism{FT},
                         PalmLikeOrganism{FT},
                         TreeLikeOrganism{FT}},
            flow::FT
) where {FT<:AbstractFloat}
    roots_flow!(plant.roots,
                plant.cache_k,
                plant.cache_p,
                plant.cache_q,
                flow);

    return nothing
end








###############################################################################
#
# Calculate root flow rates and pressure from total flow rate
#
###############################################################################
"""
    xylem_flow(root::RootHydraulics{FT},
               pressure::FT,
               ini::FT = FT(1)
    ) where {FT<:AbstractFloat}

Calculate the flow rate from a given tree base pressure, given
- `root` [`RootHydraulics`](@ref) type struct
- `pressure` Given tree base pressure in `[MPa]`
- `ini` Initial guess
"""
function xylem_flow(
            root::RootHydraulics{FT},
            pressure::FT,
            ini::FT = FT(1)
) where {FT<:AbstractFloat}
    _fh = (root.p_ups + root.p_osm * root.T_sap / T₂₅(FT) - pressure) *
          root.k_max / root.f_vis;
    _fl = -_fh;
    _fx = min(_fh, ini);
    _ms = NewtonBisectionMethod{FT}(x_min=_fl, x_max=_fh, x_ini=_fx);
    _st = SolutionTolerance{FT}(eps(FT)*100, 50);
    @inline f(x) = end_pressure(root, x) - pressure;
    _solut = find_zero(f, _ms, _st);

    return _solut
end








###############################################################################
#
# Update the flow profile from leaf flow rates
#
###############################################################################
function flow_profile!(hs::GrassLikeOrganism{FT}) where {FT<:AbstractFloat}
    # leaf rate is per leaf area so stem flow should that times leaf area
    _flow::FT = 0;
    for _i in eachindex(hs.leaves)
        _flow += hs.leaves[_i].flow * hs.leaves[_i].area;
    end

    # update root flow among root layers
    roots_flow!(hs, _flow);
    for _i in eachindex(hs.roots)
        hs.roots[_i].flow = hs.cache_q[_i];
    end

    return nothing
end




function flow_profile!(hs::PalmLikeOrganism{FT}) where {FT<:AbstractFloat}
    # leaf rate is per leaf area so stem flow should that times leaf area
    _flow::FT = 0;
    for _i in eachindex(hs.leaves)
        _flow += hs.leaves[_i].flow * hs.leaves[_i].area;
    end

    # trunk flow rate
    hs.trunk.flow = _flow;

    # update root flow among root layers
    roots_flow!(hs, _flow);
    for _i in eachindex(hs.roots)
        hs.roots[_i].flow = hs.cache_q[_i];
    end

    return nothing
end




function flow_profile!(hs::TreeLikeOrganism{FT}) where {FT<:AbstractFloat}
    # leaf rate is per leaf area so stem flow should that times leaf area
    _flow::FT = 0;
    for _i in eachindex(hs.leaves)
        hs.branch[_i].flow = hs.leaves[_i].flow * hs.leaves[_i].area;
        _flow += hs.branch[_i].flow
    end

    # trunk flow rate
    hs.trunk.flow = _flow;

    # update root flow among root layers
    roots_flow!(hs, _flow);
    for _i in eachindex(hs.roots)
        hs.roots[_i].flow = hs.cache_q[_i];
    end

    return nothing
end

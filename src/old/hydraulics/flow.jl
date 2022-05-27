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
    @inline f(x) = xylem_end_pressure(hs, x) - hs.p_crt;
    _solut  = find_zero(f, _ms, _rt);

    if isnan(_solut)
        @warn twarn("E_crit is NaN, please check the settings...") hs.p_ups;
    end

    return _solut
end




function critical_flow(
            tree::MonoElementSPAC{FT},
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
    @inline f(x) = xylem_end_pressure(tree, x) - (tree.leaf).p_crt;
    _solut  = find_zero(f, _ms, _rt);

    return _solut
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
    _fh = (root.p_ups + root.p_osm * root.T_sap / T_25(FT) - pressure) *
          root.k_max / root.f_vis;
    _fl = -_fh;
    _fx = min(_fh, ini);
    _ms = NewtonBisectionMethod{FT}(x_min=_fl, x_max=_fh, x_ini=_fx);
    _st = SolutionTolerance{FT}(eps(FT)*100, 50);
    @inline f(x) = xylem_end_pressure(root, x) - pressure;
    _solut = find_zero(f, _ms, _st);

    return _solut
end








###############################################################################
#
# Update the flow profile from leaf flow rates
#
###############################################################################
#=
function flow_profile!(hs::MonoGrassSPAC{FT}) where {FT<:AbstractFloat}
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




function flow_profile!(hs::MonoPalmSPAC{FT}) where {FT<:AbstractFloat}
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




function flow_profile!(hs::MonoTreeSPAC{FT}) where {FT<:AbstractFloat}
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
=#

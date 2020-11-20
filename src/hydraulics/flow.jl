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
- `tree` [`TreeSimple`](@refy type struct)

Note, for the safety of no NaN, update critical_flow when ΔP >= -0.01
"""
function critical_flow(
            hs::LeafHydraulics{FT},
            ini::FT = FT(0.5)
) where {FT<:AbstractFloat}
    # calculate maximal flow
    _fh     = (hs.p_ups - hs.p_crt) * hs.k_sla / hs.f_vis;
    _fl     = FT(0);
    _fx     = min((_fh+_fl)/2, ini);
    _ms     = NewtonBisectionMethod{FT}(_fl, _fh, _fx);
    _rt     = SolutionTolerance{FT}(1e-5, 50);
    @inline f(x) = end_pressure(hs, x) - hs.p_crt;
    _solut  = find_zero(f, _ms, _rt);

    if isnan(_solut)
        @warn "E_crit is NaN, please check the settings...";
        @show hs.p_ups;
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
    _ms = NewtonBisectionMethod{FT}(_fl, _fh, _fx);
    _rt = SolutionTolerance{FT}(1e-3, 50);
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
    _fh    = (root.p_ups - pressure) * root.k_max / root.f_vis;
    _fl    = -_fh;
    _fx    = min(_fh, ini);
    _ms    = NewtonBisectionMethod{FT}(_fl, _fh, _fx);
    _st    = SolutionTolerance{FT}(1e-4, 50);
    @inline f(x) = end_pressure(root, x) - pressure;
    _solut = find_zero(f, _ms, _st);

    return _solut
end

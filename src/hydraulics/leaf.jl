###############################################################################
#
# Calculate xylem end risk factor (k_ratio)
#
###############################################################################
"""
    leaf_xylem_risk(hs::LeafHydraulics{FT}, flow::FT) where {FT<:AbstractFloat}

Evaluate the hydraulic risk at the end of leaf xylem, given
- `hs` [`LeafHydraulics`](@ref) type struct
- `flow` Flow rate (per leaf area)
"""
function leaf_xylem_risk(
            hs::LeafHydraulics{FT},
            flow::FT
) where {FT<:AbstractFloat}
    @unpack f_st, f_vis, vc = hs;

    p_25 = xylem_p_from_flow(hs, flow) / hs.f_st;
    k_25 = xylem_k_ratio(vc, p_25, f_vis);

    return k_25
end








###############################################################################
#
# Calculate leaf-level e_crit
#
###############################################################################
"""
    leaf_e_crit(hs::LeafHydraulics{FT}, ini::FT) where {FT<:AbstractFloat}

Calculate the critical flow rate (K ≈ 0), given
- `hs` [`LeafHydraulics`](@ref) type struct
- `ini` Initial guess

Note, for the safety of no NaN, update e_crit when ΔP >= -0.01
"""
function leaf_e_crit(
            hs::LeafHydraulics{FT},
            ini::FT = FT(0.5)
) where {FT<:AbstractFloat}
    # calculate maximal flow
    _fh     = (hs.p_ups - hs.p_crt) * hs.k_sla / hs.f_vis;
    _fl     = FT(0);
    _fx     = min((_fh+_fl)/2, ini);
    _ms     = NewtonBisectionMethod{FT}(_fl, _fh, _fx);
    _rt     = SolutionTolerance{FT}(1e-5);
    @inline f(x) = xylem_p_from_flow(hs, x) - hs.p_crt;
    _solut  = find_zero(f, _ms, _rt);
    _ec::FT = _solut;

    if isnan(_ec)
        println("E_crit is NaN, please check the settings...")
        @show hs.p_ups;
        println("");
    end

    return _ec
end

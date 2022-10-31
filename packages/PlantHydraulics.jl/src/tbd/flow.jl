#=
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
=#

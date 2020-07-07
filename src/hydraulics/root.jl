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
    _st    = SolutionTolerance{FT}(1e-4);
    @inline f(x) = xylem_p_from_flow(root, x) - pressure;
    _solut = find_zero(f, _ms, _st);

    return _solut
end




"""
    q_diff(roots::Array{RootHydraulics{FT},1}, pressure::FT, target::FT) where {FT<:AbstractFloat}

Function to calculate the flow rate difference to find solution pressure, given
- `roots` Array of [`RootHydraulics`](@ref) type struct
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
    root_qs_p_from_q(roots::Array{RootHydraulics{FT},1}, flow::FT, ini::FT) where {FT<:AbstractFloat}

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
    _pl    = FT(-20);
    _ph    = FT(20);
    _px    = min(_pl/2, ini);
    _ms    = NewtonBisectionMethod{FT}(_pl, _ph, _px);
    _st    = ResidualTolerance{FT}(1e-5);
    @inline f(x) = q_diff(roots, x, flow);
    _solut = find_zero(f, _ms, _st);

    _qs = FT[root_q_from_pressure(root, _solut, root.flow) for root in roots];

    return _qs, _solut
end

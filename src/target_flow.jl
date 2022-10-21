#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-01: migrate the function
#     2022-Jun-01: add method for LeafHydraulics
#     2022-Jun-01: add method for MonoElementSPAC
#     2022-Oct-21: add a p_ups if statement
#
#######################################################################################################################################################################################################
"""
This function returns the critical flow rate that triggers a given amount of loss of hydraulic conductance for
- Leaf hydraulic system
- Mono element SPAC system

"""
function critical_flow end


"""

    critical_flow(hs::LeafHydraulics{FT}, T::FT, ini::FT = FT(0.5); kr::FT = FT(0.001)) where {FT<:AbstractFloat}

Return the critical flow rate that triggers a given amount of loss of conductance, given
- `hs` `LeafHydraulics` type struct
- `T` Liquid temperature
- `ini` Initial guess
- `kr` Reference conductance, default is 0.001
"""
critical_flow(hs::LeafHydraulics{FT}, T::FT, ini::FT = FT(0.5); kr::FT = FT(0.001)) where {FT<:AbstractFloat} = (
    @unpack K_SLA, VC = hs;

    # compute the misc variables
    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_crt = critical_pressure(VC, kr) * _f_st;

    # add a judgement to make sure p_ups is higher than _p_crt
    if hs.p_ups < _p_crt
        return eps(FT)
    end;

    # set up method to calculate critical flow
    _fh = (hs.p_ups - _p_crt) * K_SLA / _f_vis;
    _fl = FT(0);
    _fx = min((_fh+_fl)/2, ini);
    _ms = NewtonBisectionMethod{FT}(x_min=_fl, x_max=_fh, x_ini=_fx);
    _st = SolutionTolerance{FT}(eps(FT)*100, 50);

    # define the target function
    @inline f(x) = xylem_end_pressure(hs, x, T) - _p_crt;

    # find the solution
    _solut = find_zero(f, _ms, _st);

    # warning if the solution is NaN
    if isnan(_solut)
        @warn "E_crit is NaN, please check the settings..." hs.p_ups;
    end;

    return _solut
);


"""

    critical_flow(spac::MonoElementSPAC{FT}, ini::FT = FT(0.5); kr::FT = FT(0.001)) where {FT<:AbstractFloat}

Return the critical flow rate that triggers a given amount of loss of conductance, given
- `spac` `MonoElementSPAC` type struct
- `ini` Initial guess
- `kr` Reference conductance, default is 0.001
"""
critical_flow(spac::MonoElementSPAC{FT}, ini::FT = FT(0.5); kr::FT = FT(0.001)) where {FT<:AbstractFloat} = (
    @unpack LEAF, ROOT, STEM = spac;

    # read out the conductances
    _kr = ROOT.HS.AREA * ROOT.HS.K_X / ROOT.HS.L / relative_viscosity(ROOT.t);
    _ks = STEM.HS.AREA * STEM.HS.K_X / STEM.HS.L / relative_viscosity(STEM.t);
    _kl = LEAF.HS.K_SLA / relative_viscosity(LEAF.t) * LEAF.HS.AREA;
    _kt = 1 / (1 / _kr + 1 / _ks + 1 / _kl);

    # compute leaf critical pressure
    _p_crt = critical_pressure(LEAF.HS.VC, kr) * relative_surface_tension(LEAF.t);

    # add a judgement to make sure p_ups is higher than _p_crt
    if (ROOT.HS.p_ups < _p_crt)
        return eps(FT)
    end;

    # set up method to calculate critical flow
    _fh = -_p_crt * _kt;
    _fl = FT(0);
    _fx = min((_fh+_fl)/2, ini);
    _ms = NewtonBisectionMethod{FT}(x_min=_fl, x_max=_fh, x_ini=_fx);
    _st = SolutionTolerance{FT}(eps(FT)*100, 50);

    # define the target function
    @inline f(x) = xylem_end_pressure(spac, x) - _p_crt;

    # find the solution
    _solut  = find_zero(f, _ms, _st);

    # warning if the solution is NaN
    if isnan(_solut)
        @warn "E_crit is NaN, please check the settings..." ROOT.HS.p_ups;
    end;

    return _solut
);

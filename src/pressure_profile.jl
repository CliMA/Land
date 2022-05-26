
function xylem_end_pressure end


# flow is per leaf area
xylem_end_pressure(hs::LeafHydraulics{FT}, flow::FT, t::FT) where {FT<:AbstractFloat} = (
    @unpack K_SLA, N, VC = hs;

    _f_st = relative_surface_tension(t);
    _f_vis = relative_viscosity(t);
    _p_end::FT = hs.p_ups;

    # compute k from temperature and history, then update the pressure
    for _i in eachindex(hs.k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs.k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _k = relative_hydraulic_conductance(hs.VC, _p_25) / _f_vis * K_SLA * N;
        else
            _k = _k_mem / _f_vis * K_SLA * N;
        end;

        _p_end -= flow / _k;
    end;

    return _p_end
);



xylem_end_pressure(hs::RootHydraulics{FT}, flow::FT, t::FT) where {FT<:AbstractFloat} = (
    @unpack K_MAX, K_RHIZ, N, SH, VC, ΔH = hs;

    _f_st = relative_surface_tension(t);
    _f_vis = relative_viscosity(t);
    _p_end::FT = hs.p_ups;

    # convert pressure to that at 25 °C to compute soil water content
    _p_25 = _p_end / _f_st;

    # divide the rhizosphere component based on the conductance (each ring has the same maximum conductance)
    for _i in 1:10
        _k = relative_hydraulic_conductance(SH, true, _p_25) * K_RHIZ * 10 / _f_vis;
        _p_25 -= flow / _k;
    end;

    # convert the end pressure back to that at liquid pressure to be matric potential
    _p_end = _p_25 * _f_st + hs.ψ_osm * t / T_25(FT);

    # compute k from temperature, history, and gravity, then update pressure
    for _i in eachindex(hs.k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs.k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _k = relative_hydraulic_conductance(hs.VC, _p_25) / _f_vis * K_MAX * N;
        else
            _k = _k_mem / _f_vis * K_MAX * N;
        end;

        _p_end -= flow / _k + ρg_MPa(FT) * ΔH / N;
    end;

    return _p_end
);



xylem_end_pressure(hs::StemHydraulics{FT}, flow::FT, t::FT) where {FT<:AbstractFloat} = (
    @unpack K_MAX, N, VC, ΔH = hs;

    _f_st = relative_surface_tension(t);
    _f_vis = relative_viscosity(t);
    _p_end::FT = hs.p_ups;

    # compute k from temperature, history, and gravity, then update pressure
    for _i in eachindex(hs.k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs.k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _k = relative_hydraulic_conductance(hs.VC, _p_25) / _f_vis * K_MAX * N;
        else
            _k = _k_mem / _f_vis * K_MAX * N;
        end;

        _p_end -= flow / _k + ρg_MPa(FT) * ΔH / N;
    end;

    return _p_end
);



xylem_end_pressure(organ::Union{Leaf{FT}, Root{FT}, Stem{FT}}, flow::FT) where {FT<:AbstractFloat} = xylem_end_pressure(organ.HS, flow, organ.t);



xylem_end_pressure(spac::MonoElementSPAC{FT}, flow::FT) where {FT<:AbstractFloat} = (
    @unpack LEAF, ROOT, STEM = spac;

    # calculate the p_dos for roots
    STEM.p_ups = xylem_end_pressure(ROOT, flow);
    LEAF.p_ups = xylem_end_pressure(STEM, flow);

    return xylem_end_pressure(LEAF, flow / LEAF.HS.AREA);
);



xylem_end_pressure(spac::MonoElementSPAC{FT}, f_sl::FT, f_sh::FT, r_sl::FT) where {FT<:AbstractFloat} = (
    @unpack LEAF, ROOT, STEM = spac;

    # calculate the p_dos for roots
    STEM.p_ups = xylem_end_pressure(ROOT, f_sl + f_sh;);
    LEAF.p_ups = xylem_end_pressure(STEM, f_sl + f_sh;);

    # calculate the p_dos for leaves
    p_sl = xylem_end_pressure(LEAF, f_sl / (LEAF.HS.AREA * r_sl));
    p_sh = xylem_end_pressure(LEAF, f_sh / (LEAF.HS.AREA * (1 - r_sl)));

    return p_sl, p_sh
);

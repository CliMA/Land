
function root_pk end



root_pk(hs::RootHydraulics{FT}, mode::SteadyStateFlow{FT}, T::FT) where {FT<:AbstractFloat} = (
    @unpack K_MAX, K_RHIZ, N, SH, VC, ΔH = hs;

    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;
    _r_all::FT = 0;

    # convert pressure to that at 25 °C to compute soil water content
    _p_25 = _p_end / _f_st;

    # divide the rhizosphere component based on the conductance (each ring has the same maximum conductance)
    for _i in 1:10
        _k = relative_hydraulic_conductance(SH, true, _p_25) * K_RHIZ * 10 / _f_vis;
        _p_25 -= mode.flow / _k;
        _r_all += 1 / _k;
    end;

    # convert the end pressure back to that at liquid pressure to be matric potential
    _p_end = _p_25 * _f_st + hs.ψ_osm * T / T_25(FT);

    # compute k from temperature, history, and gravity, then update pressure
    for _i in eachindex(hs.k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs.k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _kr = relative_hydraulic_conductance(hs.VC, _p_25);
            _k = _kr / _f_vis * K_MAX * N;
        else
            _k = _k_mem / _f_vis * K_MAX * N;
        end;

        _p_end -= mode.flow / _k + ρg_MPa(FT) * ΔH / N;
        _r_all += 1 / _k;
    end;

    return _p_end, 1/_r_all
);



root_pk(hs::RootHydraulics{FT}, mode::NonSteadyStateFlow{FT}, T::FT) where {FT<:AbstractFloat} = (
    @unpack K_MAX, K_RHIZ, N, SH, VC, ΔH = hs;

    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;
    _r_all::FT = 0;

    # convert pressure to that at 25 °C to compute soil water content
    _p_25 = _p_end / _f_st;

    # divide the rhizosphere component based on the conductance (each ring has the same maximum conductance)
    for _i in 1:10
        _k = relative_hydraulic_conductance(SH, true, _p_25) * K_RHIZ * 10 / _f_vis;
        _p_25 -= mode.f_in / _k;
        _r_all += 1 / _k;
    end;

    # convert the end pressure back to that at liquid pressure to be matric potential
    _p_end = _p_25 * _f_st + hs.ψ_osm * T / T_25(FT);

    # compute k from temperature, history, and gravity, then update pressure
    for _i in eachindex(hs.k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs.k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _kr = relative_hydraulic_conductance(hs.VC, _p_25);
            _k = _kr / _f_vis * K_MAX * N;
        else
            _k = _k_mem / _f_vis * K_MAX * N;
        end;

        _p_end -= mode.f_element[_i] / _k + ρg_MPa(FT) * ΔH / N;
        _r_all += 1 / _k;
    end;

    return _p_end, 1/_r_all
);



root_pk(hs::RootHydraulics{FT}, T::FT) where {FT<:AbstractFloat} = root_pk(hs, hs.FLOW, T);










function xylem_flow_profile! end



xylem_flow_profile!(hs::Union{LeafHydraulics{FT}, RootHydraulics{FT}, StemHydraulics{FT}}, mode::SteadyStateFlow{FT}, T::FT, Δt::FT) where {FT<:AbstractFloat} = nothing;



xylem_flow_profile!(hs::LeafHydraulics{FT}, mode::NonSteadyStateFlow{FT}, T::FT, Δt::FT) where {FT<:AbstractFloat} = (
    @unpack PVC, V_MAXIMUM = hs;

    _f_vis = relative_viscosity(T);

    # compute the flow rate from capacitance buffer
    mode.f_buffer[1] = (hs.p_storage - hs.p_leaf) * capacitance_buffer(PVC) / _f_vis * V_MAXIMUM;

    # make sure the buffer rate does not drain or overflow the capacictance
    if (mode.f_buffer[1] > 0) && (hs.v_storage <= mode.f_buffer[1] * Δt)
        mode.f_buffer[1] = (hs.v_storage - eps(FT)) / Δt;
    end;

    # update storage and the tissue pressure (p_storage)
    hs.v_storage -= mode.f_buffer[1] * Δt;
    hs.p_storage = xylem_pressure(PVC, hs.v_storage/V_MAXIMUM, T);

    # update flow into the tissue
    mode.f_in = mode.f_out - mode.f_buffer[1];

    return nothing
);



xylem_flow_profile!(hs::Union{RootHydraulics{FT}, StemHydraulics{FT}}, mode::NonSteadyStateFlow{FT}, T::FT, Δt::FT) where {FT<:AbstractFloat} = (
    @unpack N, PVC, V_MAXIMUM = hs;

    _f_vis = relative_viscosity(T);

    # update storage volume and pressure per slice
    _f_sum::FT = 0;
    for _i in N:-1:1
        mode.f_buffer[_i] = (hs.p_storage[_i] - hs.p_element[_i]) * capacitance_buffer(PVC) / _f_vis * V_MAXIMUM[_i];

        # make sure the buffer rate does not drain or overflow the capacictance
        if (mode.f_buffer[_i] > 0) && (hs.v_storage[_i] <= mode.f_buffer[_i] * Δt)
            mode.f_buffer[_i] = hs.v_storage[_i] / Δt;
        end;

        mode._f_diff[_i] = _f_sum;
        hs.v_storage[_i] -= mode.f_buffer[_i] * Δt;
        hs.p_storage[_i] = xylem_pressure(PVC, hs.v_storage[_i]/V_MAXIMUM[_i], T);
        mode.f_element[_i] = mode.f_out - _f_sum;
        _f_sum += mode.f_buffer[_i];
    end;

    # update flow into the tissue
    mode.f_in = mode.f_out - _f_sum;

    return nothing
);



xylem_flow_profile!(hs::RootHydraulics{FT}, mode::SteadyStateFlow, flow::FT) where {FT<:AbstractFloat} = (
    mode.flow = flow;

    return nothing
);



xylem_flow_profile!(hs::RootHydraulics{FT}, mode::NonSteadyStateFlow, f_out::FT) where {FT<:AbstractFloat} = (
    mode.f_out = f_out;
    hs.f_element .= f_out .- mode.q_diff;

    return nothing
);



xylem_flow_profile!(hs::Union{LeafHydraulics{FT}, RootHydraulics{FT}, StemHydraulics{FT}}, T::FT, Δt::FT) where {FT<:AbstractFloat} = xylem_flow_profile!(hs, hs.FLOW, T, Δt);



xylem_flow_profile!(organ::Union{Leaf{FT}, Root{FT}, Stem{FT}}, T::FT, Δt::FT) where {FT<:AbstractFloat} = xylem_flow_profile!(organ.HS, T, Δt);



xylem_flow_profile!(roots::Vector{RootHydraulics{FT}}, cache_f::Vector{FT}, cache_k::Vector{FT}, cache_p::Vector{FT}, f_sum::FT, T::FT, Δt::FT) where {FT<:AbstractFloat} = (
    # update root buffer rates to get an initial guess (flow rate not changing now)
    xylem_flow_profile!.(roots, T, FT(0));

    # recalculate the flow profiles to make sure sum are the same as f_sum
    _count = 0;
    while _count < 20
        # sync the values to ks, ps, and qs
        for _i in eachindex(roots)
            _root = roots[_i];
            xylem_flow_profile!(roots[_i], roots[_i].FLOW, cache_f[_i]);
            cache_p[_i],cache_k[_i] = root_pk(_root, T);
        end;

        # use ps and ks to compute the Δf to adjust
        _pm = mean(cache_p);
        for _i in eachindex(roots)
            cache_f[_i] -= (_pm - cache_p[_i]) * cache_k[_i];
        end;

        # adjust the fs so that sum(fs) = f_sum
        _f_diff = sum(cache_f) - f_sum;
        if (abs(_f_diff) < FT(1e-6)) && (maximum(cache_p) - minimum(cache_p) < 1e-4)
            break
        end;
        _k_sum  = sum(cache_k);
        for _i in eachindex(roots)
            cache_f[_i] -= _f_diff * cache_k[1] / _k_sum;
        end;
    end;

    # update root buffer rates again
    for _i in eachindex(roots)
        xylem_flow_profile!(roots[_i], roots[_i].FLOW, cache_f[_i]);
        xylem_flow_profile!(roots[_i], T, Δt);
    end;

    return nothing
);



xylem_flow_profile!(spac::MonoElementSPAC{FT}) where {FT<:AbstractFloat} = (
    return nothing
);



#=
xylem_flow_profile!(spac::MonoGrassSPAC{FT}) where {FT<:AbstractFloat} = (
    # leaf rate is per leaf area so stem flow should that times leaf area
    _flow::FT = 0;
    for _i in eachindex(hs.leaves)
        _flow += hs.leaves[_i].flow * hs.leaves[_i].area;
    end;

    # update root flow among root layers
    roots_flow!(hs, _flow);
    for _i in eachindex(hs.roots)
        hs.roots[_i].flow = hs.cache_q[_i];
    end;

    return nothing
);
=#

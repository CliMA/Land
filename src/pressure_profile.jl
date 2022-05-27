#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-May-25: migrate function to version v0.3
#     2022-May-24: rename the function to xylem_end_pressure
#
#######################################################################################################################################################################################################
"""
This function returns the xylem end water pressure. The supported methods are

$(METHODLIST)

"""
function xylem_end_pressure end


#######################################################################################################################################################################################################
#
# Changes to the method
# General
#     2022-May-25: add method for LeafHydraulics
#
#######################################################################################################################################################################################################
"""

    xylem_end_pressure(hs::LeafHydraulics{FT}, flow::FT, T::FT) where {FT<:AbstractFloat}

Return the xylem end water pressure in MPa, given
- `hs` `LeafHydraulics` type struct
- `flow` Flow rate per leaf area `[mol m⁻² s⁻¹]`
- `T` Temperature
"""
xylem_end_pressure(hs::LeafHydraulics{FT}, flow::FT, T::FT) where {FT<:AbstractFloat} = (
    @unpack K_SLA, N, VC = hs;

    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
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


#######################################################################################################################################################################################################
#
# Changes to the method
# General
#     2022-May-25: add method for RootHydraulics
#
#######################################################################################################################################################################################################
"""

    xylem_end_pressure(hs::RootHydraulics{FT}, flow::FT, T::FT) where {FT<:AbstractFloat}

Return the xylem end water pressure in MPa, given
- `hs` `RootHydraulics` type struct
- `flow` Flow rate `[mol s⁻¹]`
- `T` Temperature
"""
xylem_end_pressure(hs::RootHydraulics{FT}, flow::FT, T::FT) where {FT<:AbstractFloat} = (
    @unpack K_MAX, K_RHIZ, N, SH, VC, ΔH = hs;

    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;

    # convert pressure to that at 25 °C to compute soil water content
    _p_25 = _p_end / _f_st;

    # divide the rhizosphere component based on the conductance (each ring has the same maximum conductance)
    for _i in 1:10
        _k = relative_hydraulic_conductance(SH, true, _p_25) * K_RHIZ * 10 / _f_vis;
        _p_25 -= flow / _k;
    end;

    # convert the end pressure back to that at liquid pressure to be matric potential
    _p_end = _p_25 * _f_st + hs.ψ_osm * T / T_25(FT);

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


#######################################################################################################################################################################################################
#
# Changes to the method
# General
#     2022-May-25: add method for StemHydraulics
#
#######################################################################################################################################################################################################
"""

    xylem_end_pressure(hs::StemHydraulics{FT}, flow::FT, T::FT) where {FT<:AbstractFloat}

Return the xylem end water pressure in MPa, given
- `hs` `StemHydraulics` type struct
- `flow` Flow rate `[mol s⁻¹]`
- `T` Temperature
"""
xylem_end_pressure(hs::StemHydraulics{FT}, flow::FT, T::FT) where {FT<:AbstractFloat} = (
    @unpack K_MAX, N, VC, ΔH = hs;

    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
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


#######################################################################################################################################################################################################
#
# Changes to the method
# General
#     2022-May-25: add method for Leaf, Root, and Stem
#
#######################################################################################################################################################################################################
"""

    xylem_end_pressure(organ::Union{Leaf{FT}, Root{FT}, Stem{FT}}, flow::FT) where {FT<:AbstractFloat}

Return the xylem end water pressure in MPa, given
- `organ` `Leaf`, `Root`, or `Stem` type struct
- `flow` Flow rate (per leaf area for Leaf) `[mol (m⁻²) s⁻¹]`
"""
xylem_end_pressure(organ::Union{Leaf{FT}, Root{FT}, Stem{FT}}, flow::FT) where {FT<:AbstractFloat} = xylem_end_pressure(organ.HS, flow, organ.t);


#######################################################################################################################################################################################################
#
# Changes to the method
# General
#     2022-May-25: add method for MonoElementSPAC (without partitioning)
#
#######################################################################################################################################################################################################
"""

    xylem_end_pressure(spac::MonoElementSPAC{FT}, flow::FT) where {FT<:AbstractFloat}

Return the xylem end water pressure in MPa, given
- `spac` `MonoElementSPAC` type struct
- `flow` Flow rate `[mol s⁻¹]`
"""
xylem_end_pressure(spac::MonoElementSPAC{FT}, flow::FT) where {FT<:AbstractFloat} = (
    @unpack LEAF, ROOT, STEM = spac;

    # calculate the p_dos for roots
    STEM.HS.p_ups = xylem_end_pressure(ROOT, flow);
    LEAF.HS.p_ups = xylem_end_pressure(STEM, flow);

    return xylem_end_pressure(LEAF, flow / LEAF.HS.AREA);
);


#######################################################################################################################################################################################################
#
# Changes to the method
# General
#     2022-May-25: add method for MonoElementSPAC (with sunlit and shaded partitioning)
# Todos
#     TODO: abstractize the flow rates with canopy RT module
#
#######################################################################################################################################################################################################
"""

    xylem_end_pressure(spac::MonoElementSPAC{FT}, f_sl::FT, f_sh::FT, r_sl::FT) where {FT<:AbstractFloat}

Return the xylem end water pressure in MPa for sunlit and shaded leaves, given
- `spac` `MonoElementSPAC` type struct
- `f_sl` Sunlit leaves flow rate `[mol s⁻¹]`
- `f_sh` Shaded leaves flow rate `[mol s⁻¹]`
- `r_sl` Sunlit leaves fraction
"""
xylem_end_pressure(spac::MonoElementSPAC{FT}, f_sl::FT, f_sh::FT, r_sl::FT) where {FT<:AbstractFloat} = (
    @unpack LEAF, ROOT, STEM = spac;

    # calculate the p_dos for roots
    STEM.HS.p_ups = xylem_end_pressure(ROOT, f_sl + f_sh;);
    LEAF.HS.p_ups = xylem_end_pressure(STEM, f_sl + f_sh;);

    # calculate the p_dos for leaves
    p_sl = xylem_end_pressure(LEAF, f_sl / (LEAF.HS.AREA * r_sl));
    p_sh = xylem_end_pressure(LEAF, f_sh / (LEAF.HS.AREA * (1 - r_sl)));

    return p_sl, p_sh
);



# TODO: update the flow in another function
function xylem_pressure_profile! end



xylem_pressure_profile!(hs::LeafHydraulics{FT}, mode::SteadyStateFlow{FT}, T::FT; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack K_SLA, N, VC = hs;

    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;

    # compute k from temperature and history, then update the pressure
    for _i in eachindex(hs.k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs.k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _kr = relative_hydraulic_conductance(hs.VC, _p_25);
            _k = _kr / _f_vis * K_SLA * N;
            if update
                hs.p_history[_i] = _p_25;
                hs.k_history[_i] = _kr;
            end;
        else
            _k = _k_mem / _f_vis * K_SLA * N;
        end;

        _p_end -= mode.flow / _k;

        hs.p_element[_i] = _p_end;
    end;

    # update the xylem end pressure
    hs.p_dos = _p_end;

    # TODO: update the leaf water potential based on extra-xylary conductance
    hs.p_leaf = _p_end - mode.flow / hs.K_OX;

    return nothing
);



xylem_pressure_profile!(hs::LeafHydraulics{FT}, mode::NonSteadyStateFlow{FT}, T::FT; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack K_SLA, N, VC = hs;

    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;

    # compute k from temperature and history, then update the pressure
    for _i in eachindex(hs.k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs.k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _kr = relative_hydraulic_conductance(hs.VC, _p_25);
            _k = _kr / _f_vis * K_SLA * N;
            if update
                hs.p_history[_i] = _p_25;
                hs.k_history[_i] = _kr;
            end;
        else
            _k = _k_mem / _f_vis * K_SLA * N;
        end;

        _p_end -= mode.f_in / _k;

        hs.p_element[_i] = _p_end;
    end;

    # update the xylem end pressure
    hs.p_dos = _p_end;

    # TODO: update the leaf water potential based on extra-xylary conductance
    hs.p_leaf = _p_end - mode.f_out / hs.K_OX;

    return nothing
);



xylem_pressure_profile!(hs::RootHydraulics{FT}, mode::SteadyStateFlow{FT}, T::FT; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack K_MAX, K_RHIZ, N, SH, VC, ΔH = hs;

    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;

    # convert pressure to that at 25 °C to compute soil water content
    _p_25 = _p_end / _f_st;

    # divide the rhizosphere component based on the conductance (each ring has the same maximum conductance)
    for _i in 1:10
        _k = relative_hydraulic_conductance(SH, true, _p_25) * K_RHIZ * 10 / _f_vis;
        _p_25 -= mode.flow / _k;
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
            if update
                hs.p_history[_i] = _p_25;
                hs.k_history[_i] = _kr;
            end;
        else
            _k = _k_mem / _f_vis * K_MAX * N;
        end;

        _p_end -= mode.flow / _k + ρg_MPa(FT) * ΔH / N;

        hs.p_element[_i] = _p_end;
    end;

    # update the xylem end pressure
    hs.p_dos = _p_end;

    return nothing
);



xylem_pressure_profile!(hs::RootHydraulics{FT}, mode::NonSteadyStateFlow{FT}, T::FT; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack K_MAX, K_RHIZ, N, SH, VC, ΔH = hs;

    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;

    # convert pressure to that at 25 °C to compute soil water content
    _p_25 = _p_end / _f_st;

    # divide the rhizosphere component based on the conductance (each ring has the same maximum conductance)
    for _i in 1:10
        _k = relative_hydraulic_conductance(SH, true, _p_25) * K_RHIZ * 10 / _f_vis;
        _p_25 -= mode.f_in / _k;
    end;

    # convert the end pressure back to that at liquid pressure to be matric potential
    _p_end = _p_25 * _f_st + hs.ψ_osm * T / T_25(FT);
    hs.p_rhiz = _p_end;

    # compute k from temperature, history, and gravity, then update pressure
    for _i in eachindex(hs.k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs.k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _kr = relative_hydraulic_conductance(hs.VC, _p_25);
            _k = _kr / _f_vis * K_MAX * N;
            if update
                hs.p_history[_i] = _p_25;
                hs.k_history[_i] = _kr;
            end;
        else
            _k = _k_mem / _f_vis * K_MAX * N;
        end;

        _p_end -= mode.q_element[_i] / _k + ρg_MPa(FT) * ΔH / N;

        hs.p_element[_i] = _p_end;
    end;

    # update the xylem end pressure
    hs.p_dos = _p_end;

    return nothing
);



xylem_pressure_profile!(hs::StemHydraulics{FT}, mode::SteadyStateFlow{FT}, T::FT; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack K_MAX, N, VC, ΔH = hs;

    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;

    # compute k from temperature, history, and gravity, then update pressure
    for _i in eachindex(hs.k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs.k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _kr = relative_hydraulic_conductance(hs.VC, _p_25);
            _k = _kr / _f_vis * K_MAX * N;
            if update
                hs.p_history[_i] = _p_25;
                hs.k_history[_i] = _kr;
            end;
        else
            _k = _k_mem / _f_vis * K_MAX * N;
        end;

        _p_end -= mode.flow / _k + ρg_MPa(FT) * ΔH / N;

        hs.p_element[_i] = _p_end;
    end;

    # update the xylem end pressure
    hs.p_dos = _p_end;

    return nothing
);



xylem_pressure_profile!(hs::StemHydraulics{FT}, mode::NonSteadyStateFlow{FT}, T::FT; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack K_MAX, N, VC, ΔH = hs;

    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;

    # compute k from temperature, history, and gravity, then update pressure
    for _i in eachindex(hs.k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs.k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _kr = relative_hydraulic_conductance(hs.VC, _p_25);
            _k = _kr / _f_vis * K_MAX * N;
            if update
                hs.p_history[_i] = _p_25;
                hs.k_history[_i] = _kr;
            end;
        else
            _k = _k_mem / _f_vis * K_MAX * N;
        end;

        _p_end -= mode.f_element[_i] / _k + ρg_MPa(FT) * ΔH / N;

        hs.p_element[_i] = _p_end;
    end;

    # update the xylem end pressure
    hs.p_dos = _p_end;

    return nothing
);



xylem_pressure_profile!(organ::Union{Leaf{FT}, Root{FT}, Stem{FT}}; update::Bool = true) where {FT<:AbstractFloat} = xylem_pressure_profile!(organ.HS, organ.FLOW, organ.t; update = update);



xylem_pressure_profile!(spac::MonoElementSPAC{FT}; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack LEAF, ROOT, STEM = spac;

    xylem_pressure_profile!(ROOT; update = update);
    STEM.HS.p_ups = ROOT.HS.p_dos;
    xylem_pressure_profile!(STEM; update = update);
    LEAF.HS.p_ups = STEM.HS.p_dos;
    xylem_pressure_profile!(LEAF; update = update);

    return nothing
);



xylem_pressure_profile!(spac::MonoGrassSPAC{FT}; update::Bool = false) where {FT<:AbstractFloat} = (
    @unpack LEAVES, N_ROOT, ROOTS = spac;

    # update the profile in roots
    _p_mean::FT = 0;
    for _root in ROOTS
        xylem_pressure_profile!(_root; update = update);
        _p_mean += _root.p_dos;
    end;
    _p_mean /= N_ROOT;

    # update the profile in leaves
    for _leaf in LEAVES
        _leaf.HS.p_ups = _p_mean;
        xylem_pressure_profile!(_leaf; update = update);
    end;

    return nothing
);



xylem_pressure_profile!(spac::MonoPalmSPAC{FT}; update::Bool = false) where {FT<:AbstractFloat} = (
    @unpack LEAVES, N_ROOT, ROOTS, TRUNK = spac;

    # update the profile in roots
    _p_mean::FT = 0;
    for _root in ROOTS
        xylem_pressure_profile!(_root; update = update);
        _p_mean += _root.p_dos;
    end;
    _p_mean /= N_ROOT;

    # update the profile in trunk
    TRUNK.HS.p_ups = _p_mean;
    xylem_pressure_profile!(TRUNK; update = update);

    # update the profile in leaf
    for _leaf in LEAVES
        _leaf.HS.p_ups = TRUNK.HS.p_dos;
        xylem_pressure_profile!(_leaf; update = update);
    end;

    return nothing
);



xylem_pressure_profile!(spac::MonoTreeSPAC{FT}; update::Bool = false) where {FT<:AbstractFloat} = (
    @unpack BRANCHES, LEAVES, N_ROOT, ROOTS, TRUNK = spac;

    # update the profile in roots
    _p_mean::FT = 0;
    for _root in ROOTS
        xylem_pressure_profile!(_root; update = update);
        _p_mean += _root.p_dos;
    end;
    _p_mean /= N_ROOT;

    # update the profile in trunk
    TRUNK.HS.p_ups = _p_mean;
    xylem_pressure_profile!(TRUNK; update = update);

    # update the profile in branch and leaf
    for _i in eachindex(BRANCHES)
        _stem = BRANCHES[_i];
        _leaf = LEAVES[_i];
        _stem.HS.p_ups = TRUNK.HS.p_dos;
        xylem_pressure_profile!(_stem; update = update);
        _leaf.HS.p_ups = _stem.HS.p_dos;
        xylem_pressure_profile!(_leaf; update = update);
    end;

    return nothing
);




# TODO: f_sl and f_sh are not saved correctly here, a change over Leaf structure may be required
#=
xylem_pressure_profile!(spac::MonoElementSPAC{FT}, f_sl::FT, f_sh::FT, r_sl::FT; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack LEAF, ROOT, STEM = spac;

    xylem_pressure_profile!(ROOT, f_sl + f_sh; update = update);
    STEM.HS.p_ups = ROOT.HS.p_dos;
    xylem_pressure_profile!(STEM, f_sl + f_sh; update = update);
    LEAF.HS.p_ups = STEM.HS.p_dos;
    xylem_pressure_profile!(LEAF, f_sl / (LEAF.HS.AREA * r_sl); update = update);
    xylem_pressure_profile!(LEAF, f_sh / (LEAF.HS.AREA * (1 - r_sl)); update = update);

    return nothing
);
=#

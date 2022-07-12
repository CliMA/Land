#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-May-25: migrate function to version v0.3
#     2022-May-25: rename the function to xylem_end_pressure
#     2022-May-25: add method for LeafHydraulics
#     2022-May-25: add method for RootHydraulics
#     2022-May-25: add method for StemHydraulics
#     2022-May-25: add method for Leaf, Root, and Stem
#     2022-May-25: add method for MonoElementSPAC (without partitioning)
# To do
#     TODO: add method for NonSteadyStateFlow
#
#######################################################################################################################################################################################################
"""

    xylem_end_pressure(organ::Union{Leaf{FT}, Root{FT}, Stem{FT}}, flow::FT) where {FT<:AbstractFloat}
    xylem_end_pressure(spac::MonoElementSPAC{FT}, flow::FT) where {FT<:AbstractFloat}

Return the xylem end water pressure in MPa, given
- `organ` `Leaf`, `Root`, or `Stem` type struct
- `flow` Flow rate (per leaf area for Leaf) `[mol (m⁻²) s⁻¹]`
- `spac` `MonoElementSPAC` type struct

"""
function xylem_end_pressure end

xylem_end_pressure(spac::MonoElementSPAC{FT}, flow::FT) where {FT<:AbstractFloat} = (
    @unpack LEAF, ROOT, STEM = spac;

    # calculate the p_dos for root and stem
    STEM.HS.p_ups = xylem_end_pressure(ROOT, flow);
    LEAF.HS.p_ups = xylem_end_pressure(STEM, flow);

    return xylem_end_pressure(LEAF, flow / LEAF.HS.AREA);
);

xylem_end_pressure(organ::Union{Leaf{FT}, Root{FT}, Stem{FT}}, flow::FT) where {FT<:AbstractFloat} = xylem_end_pressure(organ.HS, flow, organ.t);

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


#=
#######################################################################################################################################################################################################
#
# Changes to the method
# General
#     2022-May-25: add method for MonoElementSPAC (with sunlit and shaded partitioning)
# To do
#     TODO: abstractize the flow rates with canopy RT module
#     TODO: make this method for Leaves1D
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

    # calculate the p_dos for root and stem
    STEM.HS.p_ups = xylem_end_pressure(ROOT, f_sl + f_sh;);
    LEAF.HS.p_ups = xylem_end_pressure(STEM, f_sl + f_sh;);

    # calculate the p_dos for leaves
    p_sl = xylem_end_pressure(LEAF, f_sl / (LEAF.HS.AREA * r_sl));
    p_sh = xylem_end_pressure(LEAF, f_sh / (LEAF.HS.AREA * (1 - r_sl)));

    return p_sl, p_sh
);
=#


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-May-27: migrate function to version v0.3
#     2022-May-27: rename the function to xylem_pressure_profile!
#     2022-May-27: add method for LeafHydraulics at steady state
#     2022-May-27: add method for LeafHydraulics at non steady state
#     2022-May-27: add method for RootHydraulics at steady state
#     2022-May-27: add method for RootHydraulics at non steady state
#     2022-May-27: add method for StemHydraulics at steady state
#     2022-May-27: add method for StemHydraulics at non steady state
#     2022-May-27: add method for Leaf, Root, and Stem (for both steady and non-steady states)
#     2022-Jun-30: add support to Leaves2D
#     2022-Jun-30: add method for Leaves1D
#     2022-May-27: add method for MonoElementSPAC
#     2022-May-27: add method for MonoGrassSPAC
#     2022-May-27: add method for MonoPalmSPAC
#     2022-May-27: add method for MonoTreeSPAC
#     2022-May-31: pass the test
#     2022-Jun-29: rename SPAC to ML*SPAC to be more accurate
#     2022-Jul-08: deflate documentations
#     2022-Jul-12: compute e_crit for leaves
#     2022-Jul-12: compute β for leaves (only for empirical models)
# To do
#     TODO: add leaf extra-xylary vulnerability curve
#
#######################################################################################################################################################################################################
"""
This function is designed for the following purposes:
- Update organ pressre profile
- Update pressre profile for the entire SPAC

"""
function xylem_pressure_profile! end


"""

    xylem_pressure_profile!(organ::Union{Leaf{FT}, Leaves2D{FT}, Root{FT}, Stem{FT}}; update::Bool = true) where {FT<:AbstractFloat}
    xylem_pressure_profile!(organ::Leaves1D{FT}; update::Bool = true) where {FT<:AbstractFloat}

Update xylem pressure profile (flow profile needs to be updated a priori), given
- `organ` `Leaf`, `Leaves1D`, `Leaves2D`, `Root`, or `Stem` type organ
- `update` If true, update xylem cavitation legacy and leaf critical flow (e_crit)

"""
xylem_pressure_profile!(organ::Union{Leaf{FT}, Leaves2D{FT}, Root{FT}, Stem{FT}}; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack HS = organ;

    xylem_pressure_profile!(HS, HS.FLOW, organ.t; update = update);

    return nothing
);

xylem_pressure_profile!(organ::Leaves1D{FT}; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack HS, HS2 = organ;

    xylem_pressure_profile!(HS, HS.FLOW, organ.t[1]; update = update);
    HS2.k_history .= HS.k_history;
    HS2.p_history .= HS.p_history;

    xylem_pressure_profile!(HS2, HS2.FLOW, organ.t[2]; update = update);
    HS.k_history .= HS2.k_history;
    HS.p_history .= HS2.p_history;

    return nothing
);

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

    # update the leaf water potential based on extra-xylary conductance
    hs.p_leaf = _p_end - mode.flow / hs.K_OX;

    # update the e_crit
    hs.e_crit = critical_flow(hs, T, hs.e_crit);

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

    # update the leaf water potential based on extra-xylary conductance
    hs.p_leaf = _p_end - mode.f_out / hs.K_OX;

    # update the e_crit
    hs.e_crit = critical_flow(hs, T, hs.e_crit);

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

        _p_end -= mode.f_element[_i] / _k + ρg_MPa(FT) * ΔH / N;

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


"""

    xylem_pressure_profile!(spac::MonoElementSPAC{FT}; update::Bool = true) where {FT<:AbstractFloat}
    xylem_pressure_profile!(spac::MonoMLGrassSPAC{FT}; update::Bool = true) where {FT<:AbstractFloat}
    xylem_pressure_profile!(spac::MonoMLPalmSPAC{FT}; update::Bool = true) where {FT<:AbstractFloat}
    xylem_pressure_profile!(spac::MonoMLTreeSPAC{FT}; update::Bool = true) where {FT<:AbstractFloat}

Update xylem pressure profile (flow profile needs to be updated a priori), given
- `spac` `MonoElementSPAC`, `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` type spac
- `update` If true, update xylem cavitation legacy

"""
xylem_pressure_profile!(spac::MonoElementSPAC{FT}; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack LEAF, ROOT, STEM = spac;

    xylem_pressure_profile!(ROOT; update = update);
    STEM.HS.p_ups = ROOT.HS.p_dos;
    xylem_pressure_profile!(STEM; update = update);
    LEAF.HS.p_ups = STEM.HS.p_dos;
    xylem_pressure_profile!(LEAF; update = update);

    # update the β factor for empirical models
    β_factor!(spac);

    return nothing
);

xylem_pressure_profile!(spac::MonoMLGrassSPAC{FT}; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack LEAVES, N_ROOT, ROOTS = spac;

    # update the profile in roots
    _p_mean::FT = 0;
    for _root in ROOTS
        xylem_pressure_profile!(_root; update = update);
        _p_mean += _root.HS.p_dos;
    end;
    _p_mean /= N_ROOT;

    # update the profile in leaves
    for _leaf in LEAVES
        _leaf.HS.p_ups = _p_mean;
        xylem_pressure_profile!(_leaf; update = update);
    end;

    # update the β factor for empirical models
    β_factor!(spac);

    return nothing
);

xylem_pressure_profile!(spac::MonoMLPalmSPAC{FT}; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack LEAVES, N_ROOT, ROOTS, TRUNK = spac;

    # update the profile in roots
    _p_mean::FT = 0;
    for _root in ROOTS
        xylem_pressure_profile!(_root; update = update);
        _p_mean += _root.HS.p_dos;
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

    # update the β factor for empirical models
    β_factor!(spac);

    return nothing
);

xylem_pressure_profile!(spac::MonoMLTreeSPAC{FT}; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack BRANCHES, LEAVES, N_ROOT, ROOTS, TRUNK = spac;

    # update the profile in roots
    _p_mean::FT = 0;
    for _root in ROOTS
        xylem_pressure_profile!(_root; update = update);
        _p_mean += _root.HS.p_dos;
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

    # update the β factor for empirical models
    β_factor!(spac);

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

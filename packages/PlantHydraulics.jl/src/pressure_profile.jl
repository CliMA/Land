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
#     2022-Oct-20: use add SoilLayer to function variables, because of the removal of SH from RootHydraulics
# To do
#     TODO: add method for NonSteadyStateFlow
#
#######################################################################################################################################################################################################
"""

    xylem_end_pressure(organ::Union{Leaf{FT}, Stem{FT}}, flow::FT) where {FT<:AbstractFloat}
    xylem_end_pressure(organ::Root{FT}, slayer::SoilLayer{FT}, flow::FT) where {FT<:AbstractFloat}
    xylem_end_pressure(spac::MonoElementSPAC{FT}, flow::FT) where {FT<:AbstractFloat}

Return the xylem end water pressure in MPa, given
- `organ` `Leaf`, `Root`, or `Stem` type struct
- `slayer` Soil layer corresponded to root
- `flow` Flow rate (per leaf area for Leaf) `[mol (m⁻²) s⁻¹]`
- `spac` `MonoElementSPAC` type struct

"""
function xylem_end_pressure end

xylem_end_pressure(spac::MonoElementSPAC{FT}, flow::FT) where {FT<:AbstractFloat} = (
    @unpack LEAF, ROOT, SOIL, STEM = spac;

    # calculate the p_dos for root and stem
    STEM.HS.p_ups = xylem_end_pressure(ROOT, SOIL.LAYERS[1], flow);
    LEAF.HS.p_ups = xylem_end_pressure(STEM, flow);

    return xylem_end_pressure(LEAF, flow / LEAF.HS.AREA);
);

xylem_end_pressure(organ::Union{Leaf{FT}, Stem{FT}}, flow::FT) where {FT<:AbstractFloat} = xylem_end_pressure(organ.HS, flow, organ.t);

xylem_end_pressure(organ::Root{FT}, slayer::SoilLayer{FT}, flow::FT) where {FT<:AbstractFloat} = xylem_end_pressure(organ.HS, slayer, flow, organ.t);

xylem_end_pressure(hs::LeafHydraulics{FT}, flow::FT, T::FT) where {FT<:AbstractFloat} = (
    @unpack DIM_XYLEM, K_SLA, VC = hs;

    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;

    # compute k from temperature and history, then update the pressure
    for _i in eachindex(hs._k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs._k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _k = relative_hydraulic_conductance(hs.VC, _p_25) / _f_vis * K_SLA * DIM_XYLEM;
        else
            _k = _k_mem / _f_vis * K_SLA * DIM_XYLEM;
        end;

        _p_end -= flow / _k;
    end;

    return _p_end
);

xylem_end_pressure(hs::RootHydraulics{FT}, slayer::SoilLayer{FT}, flow::FT, T::FT) where {FT<:AbstractFloat} = (
    @unpack AREA, DIM_XYLEM, K_RHIZ, K_X, L, VC, ΔH = hs;

    _k_max = AREA * K_X / L;
    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;

    # convert pressure to that at 25 °C to compute soil water content
    _p_25 = _p_end / _f_st;

    # divide the rhizosphere component based on the conductance (each ring has the same maximum conductance)
    for _ in 1:10
        _k = relative_hydraulic_conductance(slayer.VC, true, _p_25) * K_RHIZ * 10 / _f_vis;
        _p_25 -= flow / _k;
    end;

    # convert the end pressure back to that at liquid pressure to be matric potential
    _p_end = _p_25 * _f_st + hs.ψ_osm * T / T₂₅(FT);

    # compute k from temperature, history, and gravity, then update pressure
    for _i in eachindex(hs._k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs._k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _k = relative_hydraulic_conductance(hs.VC, _p_25) / _f_vis * _k_max * DIM_XYLEM;
        else
            _k = _k_mem / _f_vis * _k_max * DIM_XYLEM;
        end;

        _p_end -= flow / _k + ρg_MPa(FT) * ΔH / DIM_XYLEM;
    end;

    return _p_end
);

xylem_end_pressure(hs::StemHydraulics{FT}, flow::FT, T::FT) where {FT<:AbstractFloat} = (
    @unpack AREA, DIM_XYLEM, K_X, L, VC, ΔH = hs;

    _k_max = AREA * K_X / L;
    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;

    # compute k from temperature, history, and gravity, then update pressure
    for _i in eachindex(hs._k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs._k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _k = relative_hydraulic_conductance(hs.VC, _p_25) / _f_vis * _k_max * DIM_XYLEM;
        else
            _k = _k_mem / _f_vis * _k_max * DIM_XYLEM;
        end;

        _p_end -= flow / _k + ρg_MPa(FT) * ΔH / DIM_XYLEM;
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
#     2022-Jul-12: compute e_crit for leaves
#     2022-Jul-12: compute β for leaves (only for empirical models)
#     2022-Jul-14: update root p_ups from SOIL
#     2022-Oct-20: use add SoilLayer to function variables, because of the removal of SH from RootHydraulics
# To do
#     TODO: add leaf extra-xylary vulnerability curve
#
#######################################################################################################################################################################################################
"""
This function is designed for the following purposes:
- Update organ pressure profile
- Update pressure profile for the entire SPAC

"""
function xylem_pressure_profile! end


"""

    xylem_pressure_profile!(organ::Union{Leaf{FT}, Leaves2D{FT}, Stem{FT}}; update::Bool = true) where {FT<:AbstractFloat}
    xylem_pressure_profile!(organ::Root{FT}, slayer::SoilLayer{FT}; update::Bool = true) where {FT<:AbstractFloat}
    xylem_pressure_profile!(organ::Leaves1D{FT}; update::Bool = true) where {FT<:AbstractFloat}

Update xylem pressure profile (flow profile needs to be updated a priori), given
- `organ` `Leaf`, `Leaves1D`, `Leaves2D`, `Root`, or `Stem` type organ
- `slayer` Soil layer corresponded to root
- `update` If true, update xylem cavitation legacy and leaf critical flow (e_crit)

"""
xylem_pressure_profile!(organ::Union{Leaf{FT}, Leaves2D{FT}, Stem{FT}}; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack HS = organ;

    xylem_pressure_profile!(HS, HS.FLOW, organ.t; update = update);

    return nothing
);

xylem_pressure_profile!(organ::Root{FT}, slayer::SoilLayer{FT}; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack HS = organ;

    xylem_pressure_profile!(HS, slayer, HS.FLOW, organ.t; update = update);

    return nothing
);

xylem_pressure_profile!(organ::Leaves1D{FT}; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack HS, HS2 = organ;

    xylem_pressure_profile!(HS, HS.FLOW, organ.t[1]; update = update);
    HS2._k_history .= HS._k_history;
    HS2.p_history .= HS.p_history;

    xylem_pressure_profile!(HS2, HS2.FLOW, organ.t[2]; update = update);
    HS._k_history .= HS2._k_history;
    HS.p_history .= HS2.p_history;

    return nothing
);

xylem_pressure_profile!(hs::LeafHydraulics{FT}, mode::SteadyStateFlow{FT}, T::FT; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack DIM_XYLEM, K_SLA, VC = hs;

    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;

    # compute k from temperature and history, then update the pressure
    for _i in eachindex(hs._k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs._k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _kr = relative_hydraulic_conductance(hs.VC, _p_25);
            _k = _kr / _f_vis * K_SLA * DIM_XYLEM;
            if update
                hs.p_history[_i] = _p_25;
                hs._k_history[_i] = _kr;
            end;
        else
            _k = _k_mem / _f_vis * K_SLA * DIM_XYLEM;
        end;

        _p_end -= mode.flow / _k;

        hs._p_element[_i] = _p_end;
    end;

    # update the xylem end pressure
    hs._p_dos = _p_end;

    # update the leaf water potential based on extra-xylary conductance
    hs.p_leaf = _p_end - mode.flow / hs.K_OX;

    # update the e_crit
    hs._e_crit = critical_flow(hs, T, hs._e_crit);

    return nothing
);

xylem_pressure_profile!(hs::LeafHydraulics{FT}, mode::NonSteadyStateFlow{FT}, T::FT; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack DIM_XYLEM, K_SLA, VC = hs;

    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;

    # compute k from temperature and history, then update the pressure
    for _i in eachindex(hs._k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs._k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _kr = relative_hydraulic_conductance(hs.VC, _p_25);
            _k = _kr / _f_vis * K_SLA * DIM_XYLEM;
            if update
                hs.p_history[_i] = _p_25;
                hs._k_history[_i] = _kr;
            end;
        else
            _k = _k_mem / _f_vis * K_SLA * DIM_XYLEM;
        end;

        _p_end -= mode.f_in / _k;

        hs._p_element[_i] = _p_end;
    end;

    # update the xylem end pressure
    hs._p_dos = _p_end;

    # update the leaf water potential based on extra-xylary conductance
    hs.p_leaf = _p_end - mode.f_out / hs.K_OX;

    # update the e_crit
    hs._e_crit = critical_flow(hs, T, hs._e_crit);

    return nothing
);

xylem_pressure_profile!(hs::RootHydraulics{FT}, slayer::SoilLayer{FT}, mode::SteadyStateFlow{FT}, T::FT; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack AREA, DIM_XYLEM, K_RHIZ, K_X, L, VC, ΔH = hs;

    _k_max = AREA * K_X / L;
    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;

    # convert pressure to that at 25 °C to compute soil water content
    _p_25 = _p_end / _f_st;

    # divide the rhizosphere component based on the conductance (each ring has the same maximum conductance)
    for _ in 1:10
        _k = relative_hydraulic_conductance(slayer.VC, true, _p_25) * K_RHIZ * 10 / _f_vis;
        _p_25 -= mode.flow / _k;
    end;

    # convert the end pressure back to that at liquid pressure to be matric potential
    _p_end = _p_25 * _f_st + hs.ψ_osm * T / T₂₅(FT);

    # compute k from temperature, history, and gravity, then update pressure
    for _i in eachindex(hs._k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs._k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _kr = relative_hydraulic_conductance(hs.VC, _p_25);
            _k = _kr / _f_vis * _k_max * DIM_XYLEM;
            if update
                hs.p_history[_i] = _p_25;
                hs._k_history[_i] = _kr;
            end;
        else
            _k = _k_mem / _f_vis * _k_max * DIM_XYLEM;
        end;

        _p_end -= mode.flow / _k + ρg_MPa(FT) * ΔH / DIM_XYLEM;

        hs._p_element[_i] = _p_end;
    end;

    # update the xylem end pressure
    hs.p_dos = _p_end;

    return nothing
);

xylem_pressure_profile!(hs::RootHydraulics{FT}, slayer::SoilLayer{FT}, mode::NonSteadyStateFlow{FT}, T::FT; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack AREA, DIM_XYLEM, K_RHIZ, K_X, L, VC, ΔH = hs;

    _k_max = AREA * K_X / L;
    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;

    # convert pressure to that at 25 °C to compute soil water content
    _p_25 = _p_end / _f_st;

    # divide the rhizosphere component based on the conductance (each ring has the same maximum conductance)
    for _ in 1:10
        _k = relative_hydraulic_conductance(slayer.VC, true, _p_25) * K_RHIZ * 10 / _f_vis;
        _p_25 -= mode.f_in / _k;
    end;

    # convert the end pressure back to that at liquid pressure to be matric potential
    _p_end = _p_25 * _f_st + hs.ψ_osm * T / T₂₅(FT);
    hs._p_rhiz = _p_end;

    # compute k from temperature, history, and gravity, then update pressure
    for _i in eachindex(hs._k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs._k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _kr = relative_hydraulic_conductance(hs.VC, _p_25);
            _k = _kr / _f_vis * _k_max * DIM_XYLEM;
            if update
                hs.p_history[_i] = _p_25;
                hs._k_history[_i] = _kr;
            end;
        else
            _k = _k_mem / _f_vis * _k_max * DIM_XYLEM;
        end;

        _p_end -= mode._f_element[_i] / _k + ρg_MPa(FT) * ΔH / DIM_XYLEM;

        hs._p_element[_i] = _p_end;
    end;

    # update the xylem end pressure
    hs.p_dos = _p_end;

    return nothing
);

xylem_pressure_profile!(hs::StemHydraulics{FT}, mode::SteadyStateFlow{FT}, T::FT; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack AREA, DIM_XYLEM, K_X, L, VC, ΔH = hs;

    _k_max = AREA * K_X / L;
    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;

    # compute k from temperature, history, and gravity, then update pressure
    for _i in eachindex(hs._k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs._k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _kr = relative_hydraulic_conductance(hs.VC, _p_25);
            _k = _kr / _f_vis * _k_max * DIM_XYLEM;
            if update
                hs.p_history[_i] = _p_25;
                hs._k_history[_i] = _kr;
            end;
        else
            _k = _k_mem / _f_vis * _k_max * DIM_XYLEM;
        end;

        _p_end -= mode.flow / _k + ρg_MPa(FT) * ΔH / DIM_XYLEM;

        hs._p_element[_i] = _p_end;
    end;

    # update the xylem end pressure
    hs.p_dos = _p_end;

    return nothing
);

xylem_pressure_profile!(hs::StemHydraulics{FT}, mode::NonSteadyStateFlow{FT}, T::FT; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack AREA, DIM_XYLEM, K_X, L, VC, ΔH = hs;

    _k_max = AREA * K_X / L;
    _f_st = relative_surface_tension(T);
    _f_vis = relative_viscosity(T);
    _p_end::FT = hs.p_ups;

    # compute k from temperature, history, and gravity, then update pressure
    for _i in eachindex(hs._k_history)
        _p_mem = hs.p_history[_i];
        _k_mem = hs._k_history[_i];

        _p_25 = _p_end / _f_st;
        if _p_25 < _p_mem
            _kr = relative_hydraulic_conductance(hs.VC, _p_25);
            _k = _kr / _f_vis * _k_max * DIM_XYLEM;
            if update
                hs.p_history[_i] = _p_25;
                hs._k_history[_i] = _kr;
            end;
        else
            _k = _k_mem / _f_vis * _k_max * DIM_XYLEM;
        end;

        _p_end -= mode._f_element[_i] / _k + ρg_MPa(FT) * ΔH / DIM_XYLEM;

        hs._p_element[_i] = _p_end;
    end;

    # update the xylem end pressure
    hs.p_dos = _p_end;

    return nothing
);


# TODO: make sure to not mixing with top soil that is meant for evaporation
# TODO: add soil ion concentration
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
    @unpack LEAF, ROOT, SOIL, STEM = spac;

    # update water potential from SOIL
    ROOT.HS.p_ups = soil_ψ_25(SOIL.LAYERS[1].VC, SOIL.LAYERS[1].θ) * relative_surface_tension(SOIL.LAYERS[1].t);

    # update pressure profiles for organs
    xylem_pressure_profile!(ROOT, SOIL.LAYERS[1]; update = update);
    STEM.HS.p_ups = ROOT.HS.p_dos;
    xylem_pressure_profile!(STEM; update = update);
    LEAF.HS.p_ups = STEM.HS.p_dos;
    xylem_pressure_profile!(LEAF; update = update);

    # update the β factor for empirical models
    β_factor!(spac);

    return nothing
);

xylem_pressure_profile!(spac::MonoMLGrassSPAC{FT}; update::Bool = true) where {FT<:AbstractFloat} = (
    @unpack DIM_ROOT, LEAVES, ROOTS, ROOTS_INDEX, SOIL = spac;

    # update water potential from SOIL
    for _i in eachindex(ROOTS_INDEX)
        ROOTS[_i].HS.p_ups = soil_ψ_25(SOIL.LAYERS[ROOTS_INDEX[_i]].VC, SOIL.LAYERS[ROOTS_INDEX[_i]].θ) * relative_surface_tension(SOIL.LAYERS[ROOTS_INDEX[_i]].t);
    end;

    # update the profile in roots
    _p_mean::FT = 0;
    for _i in eachindex(ROOTS_INDEX)
        _root = ROOTS[_i];
        _slayer = SOIL.LAYERS[ROOTS_INDEX[_i]];
        xylem_pressure_profile!(_root, _slayer; update = update);
        _p_mean += _root.HS.p_dos;
    end;
    _p_mean /= DIM_ROOT;

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
    @unpack DIM_ROOT, LEAVES, ROOTS, ROOTS_INDEX, SOIL, TRUNK = spac;

    # update water potential from SOIL
    for _i in eachindex(ROOTS_INDEX)
        ROOTS[_i].HS.p_ups = soil_ψ_25(SOIL.LAYERS[ROOTS_INDEX[_i]].VC, SOIL.LAYERS[ROOTS_INDEX[_i]].θ) * relative_surface_tension(SOIL.LAYERS[ROOTS_INDEX[_i]].t);
    end;

    # update the profile in roots
    _p_mean::FT = 0;
    for _i in eachindex(ROOTS_INDEX)
        _root = ROOTS[_i];
        _slayer = SOIL.LAYERS[ROOTS_INDEX[_i]];
        xylem_pressure_profile!(_root, _slayer; update = update);
        _p_mean += _root.HS.p_dos;
    end;
    _p_mean /= DIM_ROOT;

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
    @unpack BRANCHES, DIM_ROOT, LEAVES, ROOTS, ROOTS_INDEX, SOIL, TRUNK = spac;

    # update water potential from SOIL
    for _i in eachindex(ROOTS_INDEX)
        ROOTS[_i].HS.p_ups = soil_ψ_25(SOIL.LAYERS[ROOTS_INDEX[_i]].VC, SOIL.LAYERS[ROOTS_INDEX[_i]].θ) * relative_surface_tension(SOIL.LAYERS[ROOTS_INDEX[_i]].t);
    end;

    # update the profile in roots
    _p_mean::FT = 0;
    for _i in eachindex(ROOTS_INDEX)
        _root = ROOTS[_i];
        _slayer = SOIL.LAYERS[ROOTS_INDEX[_i]];
        xylem_pressure_profile!(_root, _slayer; update = update);
        _p_mean += _root.HS.p_dos;
    end;
    _p_mean /= DIM_ROOT;

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

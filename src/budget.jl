#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-13: add function for water budget
#     2022-Jun-14: use K_MAX and ΔZ and remove K_REF
#     2022-Jun-14: rescale rain for layer 1
#     2022-Jun-14: use METEO.rain
#     2022-Jun-14: add function for soil energy budget
#     2022-Jun-14: use METEO.rain and METEO.t_precip
#     2022-Jun-14: add net radiation energy to top soil
#     2022-Jun-15: add controller to make sure soil layers do not over saturate
#     2022-Jun-15: merge the soil_water! and soil_energy! to soil_budget!
#
#######################################################################################################################################################################################################
"""
This function have two major features:
- Compute the marginal change of soil water content and energy
- Update soil water content and energy without over-saturating or draining the soil

"""
function soil_budget! end


"""

    soil_budget!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat}

Update the marginal increase of soil water content and energy per layer, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` SPAC

"""
soil_budget!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat} = (
    @unpack METEO, ROOTS, ROOTS_INDEX, SOIL = spac;
    LAYERS = SOIL.LAYERS;

    # update soil k, ψ, and soil cp (per m³), λ_thermal for each soil layer
    for _i in 1:SOIL.N_LAYER
        LAYERS[_i].k         = relative_hydraulic_conductance(LAYERS[_i].VC, LAYERS[_i].θ) * LAYERS[_i].K_MAX * relative_viscosity(LAYERS[_i].t) / LAYERS[_i].ΔZ;
        LAYERS[_i].ψ         = soil_ψ_25.(LAYERS[_i].VC, LAYERS[_i].θ) .* relative_surface_tension.(LAYERS[_i].t);
        LAYERS[_i].cp        = LAYERS[_i].CP * LAYERS[_i].Ρ + LAYERS[_i].θ * CP_L(FT) * ρ_H₂O(FT);
        LAYERS[_i].λ_thermal = (LAYERS[_i].Λ_THERMAL + LAYERS[_i].θ * Λ_THERMAL_H₂O(FT)) / LAYERS[_i].ΔZ;
        LAYERS[_i].∂e∂t = 0;
        LAYERS[_i].∂θ∂t = 0;
    end;

    # update k, δψ, and flow rate among layers
    LAYERS[1].∂θ∂t += METEO.rain / LAYERS[1].ΔZ;
    LAYERS[1].∂e∂t += METEO.rain / LAYERS[1].ΔZ * METEO.t_precip;
    LAYERS[1].∂e∂t += SOIL.ALBEDO.r_net_lw + SOIL.ALBEDO.r_net_sw;
    for _i in 1:SOIL.N_LAYER-1
        SOIL._k[_i]         = 1 / (2 / LAYERS[_i].k + 2 / LAYERS[_i+1].k);
        SOIL._δψ[_i]        = LAYERS[_i].ψ - LAYERS[_i+1].ψ + ρg_MPa(FT) * (LAYERS[_i].Z - LAYERS[_i+1].Z);
        SOIL._q[_i]         = SOIL._k[_i] * SOIL._δψ[_i];
        SOIL._λ_thermal[_i] = 1 / (2 / LAYERS[_i].λ_thermal + 2 / LAYERS[_i+1].λ_thermal);
        SOIL._δt[_i]        = LAYERS[_i].t - LAYERS[_i+1].t;
        SOIL._q_thermal[_i] = SOIL._λ_thermal[_i] * SOIL._δt[_i];

        # if flow into the lower > 0, but the lower layer is already saturated, set the flow to 0
        if (SOIL._q[_i] > 0) && (SOIL.LAYERS[_i+1].θ >= SOIL.LAYERS[_i+1].VC.Θ_MAX)
            SOIL._q[_i] = 0;
        end;

        # if flow into the lower < 0, but the upper layer is already saturated, set the flow to 0
        if (SOIL._q[_i] < 0) && (SOIL.LAYERS[_i].θ >= SOIL.LAYERS[_i].VC.Θ_MAX)
            SOIL._q[_i] = 0;
        end;

        LAYERS[_i  ].∂θ∂t -= SOIL._q[_i] / LAYERS[_i].ΔZ;
        LAYERS[_i+1].∂θ∂t += SOIL._q[_i] / LAYERS[_i+1].ΔZ;
        LAYERS[_i  ].∂e∂t -= SOIL._q_thermal[_i] / LAYERS[_i].ΔZ;
        LAYERS[_i+1].∂e∂t += SOIL._q_thermal[_i] / LAYERS[_i+1].ΔZ;
        LAYERS[_i  ].∂e∂t -= SOIL._q[_i] * LAYERS[_i].t / LAYERS[_i].ΔZ;
        LAYERS[_i+1].∂e∂t += SOIL._q[_i] * LAYERS[_i].t / LAYERS[_i+1].ΔZ;
    end;

    # loop through the roots and compute the source/sink terms
    for _i in eachindex(ROOTS)
        LAYERS[ROOTS_INDEX[_i]].∂θ∂t -= root_sink(ROOTS[_i]) / SOIL.AREA / LAYERS[ROOTS_INDEX[_i]].ΔZ;
        LAYERS[ROOTS_INDEX[_i]].∂e∂t -= root_sink(ROOTS[_i]) / SOIL.AREA / LAYERS[ROOTS_INDEX[_i]].ΔZ * LAYERS[_i].t;
    end;

    return nothing
);


"""

    soil_budget!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, δt::FT) where {FT<:AbstractFloat}

Run soil water and energy budget, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` SPAC
- `δt` Time step

"""
soil_budget!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, δt::FT) where {FT<:AbstractFloat} = (
    @unpack SOIL = spac;

    # set runoff to 0
    SOIL.runoff = 0;

    # run the update function until time elapses
    _t_res = δt;
    while true
        _δt = adjusted_time(spac, _t_res);
        for _i in 1:SOIL.N_LAYER
            SOIL.LAYERS[_i].θ += SOIL.LAYERS[_i].∂θ∂t * _δt;
            SOIL.LAYERS[_i].e += SOIL.LAYERS[_i].∂e∂t * _δt;
        end;

        # compute surface runoff
        if SOIL.LAYERS[1].θ > SOIL.LAYERS[1].VC.Θ_SAT
            SOIL.runoff += (SOIL.LAYERS[1].θ - SOIL.LAYERS[1].VC.Θ_SAT) * SOIL.LAYERS[1].ΔZ;
            SOIL.LAYERS[1].θ = SOIL.LAYERS[1].VC.Θ_SAT;
        end;

        _t_res -= _δt;

        # if _t_res > 0 rerun the budget function, else break
        _t_res > 0 ? soil_budget!(spac) : break;
    end;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-15: add function to make sure the soil layers do not over saturate or drain
#
#######################################################################################################################################################################################################
"""

    adjusted_time(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, δt::FT) where {FT<:AbstractFloat}

Return adjusted time that soil does not over saturate or drain, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` SPAC
- `δt` Time step

"""
function adjusted_time(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, δt::FT) where {FT<:AbstractFloat}
    @unpack SOIL = spac;

    # make sure each layer does not over-saturate or drain
    _δt = δt;
    for _i in 1:SOIL.N_LAYER
        # if top soil is saturated and there is rain, _δt will not change (the rain will be counted as runoff)
        if (SOIL.LAYERS[_i].∂θ∂t > 0) && (SOIL.LAYERS[_i].θ < SOIL.LAYERS[_i].VC.Θ_SAT)
            _δt_sat = (SOIL.LAYERS[_i].VC.Θ_SAT - SOIL.LAYERS[_i].θ) / SOIL.LAYERS[_i].∂θ∂t;
            _δt = min(_δt_sat, _δt);
        elseif SOIL.LAYERS[_i].∂θ∂t < 0
            _δt_dra = (SOIL.LAYERS[_i].VC.Θ_RES - SOIL.LAYERS[_i].θ) / SOIL.LAYERS[_i].∂θ∂t;
            _δt = min(_δt_dra, _δt);
        end;
    end;

    return _δt
end


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-13: add a utility function to read root water sink
#
#######################################################################################################################################################################################################
"""

    root_sink(root::Root{FT}) where {FT<:AbstractFloat}

Return root water update, given
- `root` `Root` type struct that may contain non- and steady state flow

"""
function root_sink end

root_sink(root::Root{FT}) where {FT<:AbstractFloat} = root_sink(root.HS.FLOW);

root_sink(mode::SteadyStateFlow{FT}) where {FT<:AbstractFloat} = mode.flow;

root_sink(mode::NonSteadyStateFlow{FT}) where {FT<:AbstractFloat} = mode.f_in;

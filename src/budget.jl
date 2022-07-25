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
#     2022-Jun-16: move time stepper controller to SoilPlantAirContinuum.jl
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

    # update soil k, ψ, and λ_thermal for each soil layer
    for _i in 1:SOIL.DIM_SOIL
        LAYERS[_i].k          = relative_hydraulic_conductance(LAYERS[_i].VC, LAYERS[_i].θ) * LAYERS[_i].K_MAX * relative_viscosity(LAYERS[_i].t) / LAYERS[_i].ΔZ;
        LAYERS[_i].ψ          = soil_ψ_25(LAYERS[_i].VC, LAYERS[_i].θ) * relative_surface_tension(LAYERS[_i].t);
        LAYERS[_i]._λ_thermal = (LAYERS[_i].Λ_THERMAL + LAYERS[_i].θ * Λ_THERMAL_H₂O(FT)) / LAYERS[_i].ΔZ;
        LAYERS[_i].∂e∂t = 0;
        LAYERS[_i].∂θ∂t = 0;
    end;

    # update k, δψ, and flow rate among layers
    LAYERS[1].∂θ∂t += METEO.rain / LAYERS[1].ΔZ;
    LAYERS[1].∂e∂t += METEO.rain / LAYERS[1].ΔZ * METEO.t_precip;
    LAYERS[1].∂e∂t += SOIL.ALBEDO.r_net_lw + SOIL.ALBEDO.r_net_sw;
    for _i in 1:SOIL.DIM_SOIL-1
        SOIL._k[_i]         = 1 / (2 / LAYERS[_i].k + 2 / LAYERS[_i+1].k);
        SOIL._δψ[_i]        = LAYERS[_i].ψ - LAYERS[_i+1].ψ + ρg_MPa(FT) * (LAYERS[_i].Z - LAYERS[_i+1].Z);
        SOIL._q[_i]         = SOIL._k[_i] * SOIL._δψ[_i];
        SOIL._λ_thermal[_i] = 1 / (2 / LAYERS[_i]._λ_thermal + 2 / LAYERS[_i+1]._λ_thermal);
        SOIL._δt[_i]        = LAYERS[_i].t - LAYERS[_i+1].t;
        SOIL._q_thermal[_i] = SOIL._λ_thermal[_i] * SOIL._δt[_i];

        # if flow into the lower > 0, but the lower layer is already saturated, set the flow to 0
        if (SOIL._q[_i] > 0) && (SOIL.LAYERS[_i+1].θ >= SOIL.LAYERS[_i+1].VC.Θ_SAT)
            SOIL._q[_i] = 0;
        end;

        # if flow into the lower < 0, but the upper layer is already saturated, set the flow to 0
        if (SOIL._q[_i] < 0) && (SOIL.LAYERS[_i].θ >= SOIL.LAYERS[_i].VC.Θ_SAT)
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

    # run the time step
    for _i in 1:SOIL.DIM_SOIL
        SOIL.LAYERS[_i].θ += SOIL.LAYERS[_i].∂θ∂t * δt;
        SOIL.LAYERS[_i].e += SOIL.LAYERS[_i].∂e∂t * δt / SOIL.LAYERS[1].ΔZ;
    end;

    # compute surface runoff
    SOIL.runoff = 0;
    if SOIL.LAYERS[1].θ > SOIL.LAYERS[1].VC.Θ_SAT
        # compute top soil temperature and top soil energy out due to runoff
        _cp = SOIL.LAYERS[1].CP * SOIL.LAYERS[1].ρ + SOIL.LAYERS[1].θ * ρ_H₂O(FT) * CP_L(FT) + SOIL.runoff / SOIL.LAYERS[1].ΔZ * M_H₂O(FT) * CP_L(FT);
        _t  = SOIL.LAYERS[1].e / _cp;
        SOIL.runoff = (SOIL.LAYERS[1].θ - SOIL.LAYERS[1].VC.Θ_SAT) * SOIL.LAYERS[1].ΔZ * ρ_H₂O(FT) / M_H₂O(FT);
        SOIL.LAYERS[1].θ = SOIL.LAYERS[1].VC.Θ_SAT;
        SOIL.LAYERS[1].e -= SOIL.runoff / SOIL.LAYERS[1].ΔZ * M_H₂O(FT) * CP_L(FT) * _t;
    end;

    # update soil temperature at each layer (top layer t will be same as _t above)
    for _i in 1:SOIL.DIM_SOIL
        SOIL.LAYERS[_i]._cp = SOIL.LAYERS[_i].CP * SOIL.LAYERS[_i].ρ + SOIL.LAYERS[_i].θ * ρ_H₂O(FT) * CP_L(FT);
        SOIL.LAYERS[_i].t  = SOIL.LAYERS[_i].e / SOIL.LAYERS[_i]._cp;
    end;

    return nothing
);


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

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-13: add function for water budget
#     2022-Jun-14: use K_MAX and ΔZ and remove K_REF
#     2022-Jun-14: rescale rain for layer 1
#
#######################################################################################################################################################################################################
"""

    soil_water!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, rain::FT) where {FT<:AbstractFloat}

Update the marginal increase of soil water content per layer, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` SPAC
- `rain` Water input from top soil `[mol m⁻² s⁻¹]`

"""
function soil_water! end

soil_water!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, rain::FT) where {FT<:AbstractFloat} = (
    @unpack ROOTS, ROOTS_INDEX, SOIL = spac;
    LAYERS = SOIL.LAYERS;

    # update soil k and ψ for each soil layer
    for _i in 1:SOIL.N_LAYER
        LAYERS[_i].k = relative_hydraulic_conductance(LAYERS[_i].VC, LAYERS[_i].θ) * LAYERS[_i].K_MAX * relative_viscosity(LAYERS[_i].t) / LAYERS[_i].ΔZ;
        LAYERS[_i].ψ = soil_ψ_25.(LAYERS[_i].VC, LAYERS[_i].θ) .* relative_surface_tension.(LAYERS[_i].t);
        LAYERS[_i].∂θ∂t = 0;
    end;

    # update k, δψ, and flow rate among layers
    LAYERS[1].∂θ∂t = rain / LAYERS[1].ΔZ;
    for _i in 1:SOIL.N_LAYER-1
        SOIL._k[_i] = 1 / (2 / LAYERS[_i].k + 2 / LAYERS[_i+1].k);
        SOIL._δψ[_i] = LAYERS[_i].ψ - LAYERS[_i+1].ψ + ρg_MPa(FT) * (LAYERS[_i].Z - LAYERS[_i+1].Z);
        SOIL._q[_i] = SOIL._k[_i] * SOIL._δψ[_i];
        LAYERS[_i  ].∂θ∂t -= SOIL._q[_i] / LAYERS[_i].ΔZ;
        LAYERS[_i+1].∂θ∂t += SOIL._q[_i] / LAYERS[_i+1].ΔZ;
    end;

    # loop through the roots and compute the source/sink terms
    for _i in eachindex(ROOTS)
        LAYERS[ROOTS_INDEX[_i]].∂θ∂t -= root_sink(ROOTS[_i]) / SOIL.AREA / LAYERS[ROOTS_INDEX[_i]].ΔZ;
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


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-14: add function for soil energy budget
#
#######################################################################################################################################################################################################
"""

    soil_energy!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, rain::FT, t_rain::FT) where {FT<:AbstractFloat}

Update the marginal increase of soil energy per layer, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` SPAC
- `rain` Water input from top soil `[mol m⁻² s⁻¹]`
- `t_rain` Water temperature from top soil `[K]`

"""
function soil_energy! end

soil_energy!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, rain::FT, t_rain::FT) where {FT<:AbstractFloat} = (
    @unpack ROOTS, ROOTS_INDEX, SOIL = spac;
    LAYERS = SOIL.LAYERS;

    # update soil cp (per m³) and λ_thermal for each soil layer
    for _i in SOIL.N_LAYER
        LAYERS[_i].cp = LAYERS[_i].CP * LAYERS[_i].Ρ + LAYERS[_i].θ * CP_L(FT) * ρ_H₂O(FT);
        LAYERS[_i].λ_thermal = (LAYERS[_i].Λ_THERMAL + LAYERS[_i].θ * Λ_THERMAL_H₂O(FT)) / LAYERS[_i].ΔZ;
        LAYERS[_i].∂e∂t = 0;
    end;

    # update k, δt, and flow rate among layers (upper - lower)
    LAYERS[1].∂e∂t = rain / LAYERS[1].ΔZ * t_rain;
    for _i in SOIL.N_LAYER-1
        SOIL._λ_thermal[_i] = 1 / (2 / LAYERS[_i].λ_thermal + 2 / LAYERS[_i+1].λ_thermal);
        SOIL._δt[_i] = LAYERS[_i].t - LAYERS[_i+1].t;
        SOIL._q_thermal[_i] = SOIL._λ_thermal[_i] * SOIL._δt[_i];
        LAYERS[_i  ].∂e∂t -= SOIL._q_thermal[_i] / LAYERS[_i].ΔZ;
        LAYERS[_i+1].∂e∂t += SOIL._q_thermal[_i] / LAYERS[_i+1].ΔZ;
        LAYERS[_i  ].∂e∂t -= SOIL._q[_i] * LAYERS[_i].t / LAYERS[_i].ΔZ;
        LAYERS[_i+1].∂e∂t += SOIL._q[_i] * LAYERS[_i].t / LAYERS[_i+1].ΔZ;
    end;

    # loop through the roots and compute the source/sink terms
    for _i in eachindex(ROOTS)
        LAYERS[ROOTS_INDEX[_i]].∂e∂t -= root_sink(ROOTS[_i]) / SOIL.AREA / LAYERS[ROOTS_INDEX[_i]].ΔZ * LAYERS[_i].t;
    end;

    return nothing
);

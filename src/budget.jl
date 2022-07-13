#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-13: add a utility function to read root water sink
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

    # update soil k and ψ for each soil layer
    for _i in 1:SOIL.N_LAYER
        SOIL.LAYERS[_i].k = relative_hydraulic_conductance(SOIL.LAYERS[_i].VC, SOIL.LAYERS[_i].θ) * SOIL.LAYERS[_i].K_REF * relative_viscosity(SOIL.LAYERS[_i].t);
        SOIL.LAYERS[_i].ψ = soil_ψ_25.(SOIL.LAYERS[_i].VC, SOIL.LAYERS[_i].θ) .* relative_surface_tension.(SOIL.LAYERS[_i].t);
        SOIL.LAYERS[_i].∂θ∂t = 0;
    end;

    # update k, δψ, and flow rate among layers
    SOIL.LAYERS[1].∂θ∂t = rain;
    for _i in 1:SOIL.N_LAYER-1
        SOIL._k[_i] = 1 / (2 / SOIL.LAYERS[_i].k + 2 / SOIL.LAYERS[_i+1].k);
        SOIL._δψ[_i] = SOIL.LAYERS[_i].ψ - SOIL.LAYERS[_i+1].ψ + ρg_MPa(FT) * (SOIL.LAYERS[_i].Z - SOIL.LAYERS[_i+1].Z);
        SOIL._q[_i] = SOIL._k[_i] * SOIL._δψ[_i];
        SOIL.LAYERS[_i].∂θ∂t -= SOIL._q[_i] / (SOIL.LAYERS[_i].ZS[1] - SOIL.LAYERS[_i].ZS[2]);
        SOIL.LAYERS[_i+1].∂θ∂t += SOIL._q[_i] / (SOIL.LAYERS[_i+1].ZS[1] - SOIL.LAYERS[_i+1].ZS[2]);
    end;

    # loop through the roots and compute the source/sink terms
    for _i in eachindex(ROOTS)
        SOIL.LAYERS[ROOTS_INDEX[_i]].∂θ∂t -= root_sink(ROOTS[_i]) / SOIL.AREA / (SOIL.LAYERS[_i+1].ZS[1] - SOIL.LAYERS[_i+1].ZS[2]);
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

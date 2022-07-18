#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-15: add function to make sure the soil layers do not over saturate or drain
#     2022-Jun-18: move function from SoilHydraulics.jl to SoilPlantAirContinuum.jl
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
#     2022-Jun-18: move function from SoilHydraulics.jl to SoilPlantAirContinuum.jl
#     2022-Jun-18: make it a separate function
#
#######################################################################################################################################################################################################
"""

    time_stepper!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, δt::FT; update::Bool = false) where {FT<:AbstractFloat}

Move forward in time for SPAC with time stepper controller, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` SPAC
- `δt` Time step
- `update` If true, update leaf xylem legacy effect

"""
function time_stepper!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, δt::FT; update::Bool = false) where {FT<:AbstractFloat}
    @unpack CANOPY, LEAVES, RAD_LW, SOIL = spac;

    # run the update function until time elapses
    _t_res = δt;
    while true
        _δt = adjusted_time(spac, _t_res);

        # run the budgets for all ∂x∂t
        soil_budget!(spac, _δt);
        stomatal_conductance!(spac, _δt);
        plant_energy!(spac, _δt);
        xylem_flow_profile!(spac, _δt);

        _t_res -= _δt;

        # if _t_res > 0 rerun the budget functions (shortwave radiation not included) and etc., else break
        if _t_res > 0
            canopy_radiation!(CANOPY, LEAVES, RAD_LW, SOIL);
            xylem_pressure_profile!(spac; update = update);
            leaf_photosynthesis!(spac, GCO₂Mode());
            soil_budget!(spac);
            stomatal_conductance!(spac);
            plant_energy!(spac);
        else
            break;
        end;
    end;

    return nothing
end

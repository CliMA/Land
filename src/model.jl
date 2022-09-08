#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-12: add function to run soil plant air continuum at a time step
#     2022-Jul-13: use soil_energy! and soil_water! for soil water and energy budget
#     2022-Aug-11: run fluorescence model after photosynthesis
#     2022-Aug-18: add option θ_on to enable/disable soil water budget
#     2022-Sep-07: add method to solve for steady state solution
#
#######################################################################################################################################################################################################
"""
This function runs the model using the following steps:
- Run canopy RT model
- Run hydraulic model
- Run photosynthesis model
- Run canopy fluorescence model
- Run soil water and energy budget (calculate ∂Θ∂t and ∂e∂t only)
- Run leaf stomatal conductances (calculate ∂g∂t only)
- Run leaf energy budget (calculate ∂T∂t only)
- Run time stepper (using ∂X∂t * δt, and make sure δt is not too high)

This function is supposed to have the highest hierarchy, and should support all SPAC types defined in ClimaCache.jl. Note to update the water flow profile when initializing the SPAC.

"""
function soil_plant_air_continuum! end


# TODO: add lite mode later to update energy balance (only longwave radiation and soil+leaf energy budgets)? Or use shorter time steps (will be time consuming, but more accurate)
# TODO: add top soil evaporation
"""

    soil_plant_air_continuum!(spac::Union{MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC{FT}}, δt::FT; update::Bool = false, θ_on::Bool = true) where {FT<:AbstractFloat}
    soil_plant_air_continuum!(spac::Union{MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC{FT}}; update::Bool = false) where {FT<:AbstractFloat}

Run SPAC model and move forward in time with time stepper controller, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` SPAC
- `δt` Time step (if not given, solve for steady state solution)
- `update` If true, update leaf xylem legacy effect
- `θ_on` If true, soil water budget is on (set false to run sensitivity analysis)

"""
soil_plant_air_continuum!(spac::Union{MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC{FT}}, δt::FT; update::Bool = false, θ_on::Bool = true) where {FT<:AbstractFloat} = (
    # 1. run canopy RT
    canopy_radiation!(spac);

    # 2. run plant hydraulic model (must be run before leaf_photosynthesis! as the latter may need β for empirical models)
    xylem_pressure_profile!(spac; update = update);

    # 3. run photosynthesis model
    leaf_photosynthesis!(spac, GCO₂Mode());

    # 4. run canopy fluorescence
    canopy_fluorescence!(spac);

    # save the result at this stage for the results at the beginning of this time step

    # 5. run soil energy water budget
    soil_budget!(spac);

    # 6. run leaf stomatal conductance budget
    stomatal_conductance!(spac);

    # 7. run plant energy budget
    plant_energy!(spac);

    # 8. update the prognostic variables
    time_stepper!(spac, δt; update = update, θ_on = θ_on);

    return nothing
);

soil_plant_air_continuum!(spac::Union{MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC{FT}}; update::Bool = false) where {FT<:AbstractFloat} = (
    # 1. run canopy RT
    canopy_radiation!(spac);

    # 2. update the prognostic variables (except for soil water and temperature)
    time_stepper!(spac; update = update);

    # save the result at this stage for the results at the steady state

    return nothing
);

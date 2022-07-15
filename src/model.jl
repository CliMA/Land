#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-12: add function to run soil plant air continuum at a time step
#     2022-Jul-13: use soil_energy! and soil_water! for soil water and energy budget
#
#######################################################################################################################################################################################################
"""
This function runs the model using the following steps:
- Run canopy RT model
- Run hydraulic model
- Run photosynthesis model
- Run soil water and ebergy budget (calculate ∂Θ∂t and ∂e∂t only)
- Run leaf stomatal conductances (calculate ∂g∂t only)
- Run leaf energy budget (calculate ∂T∂t only)
- Update the prognostic variables (using ∂X∂t * δt)
- Update xylem flow profiles

This function is supposed to have the highest hierarchy, and should support all SPAC types defined in ClimaCache.jl. Note to update the water flow profile when initializing the SPAC.

"""
function soil_plant_air_continuum! end

# TODO: add lite mode later to update energy balance (only longwave radiation and soil+leaf energy budgets)? Or use shorter time steps (will be time consuming, but more accurate)
soil_plant_air_continuum!(spac::MonoMLTreeSPAC{FT}, δt::FT; update::Bool = false) where {FT<:AbstractFloat} = (
    # 1. run canopy RT
    canopy_radiation!(spac);

    # 2. run plant hydraulic model
    xylem_pressure_profile!(spac; update = update);

    # 3. run photosynthesis model
    leaf_photosynthesis!(spac, GCO₂Mode());

    # 4. run soil water budget (TODO: add top soil evaporation)
    soil_budget!(spac);

    # 5. run leaf stomatal conductance
    stomatal_conductance!(spac);

    # 6. run plant energy budget
    plant_energy!(spac);

    # 7. update the prognostic variables (TODO: update flow temperature in TRUNK, BRANCK, and LEAVES)
    soil_budget!(spac, δt);
    stomatal_conductance!(spac, δt);
    # update leaf temperature

    # 8. update xylem flow profiles from stomatal conductance
    xylem_flow_profile!(spac, δt);

    return nothing
);

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
- Run soil water budget (calculate ∂Θ∂t only)
- Run soil energy budget (calculate ∂T∂t only)
- Run leaf stomatal conductances (calculate ∂g∂t only)
- Run leaf energy budget (calculate ∂T∂t only)
- Update the prognostic variables (using ∂X∂t * δt)
- Update xylem flow profiles

This function is supposed to have the highest hierarchy, and should support all SPAC types defined in ClimaCache.jl. Note to update the water flow profile when initializing the SPAC.

"""
function soil_plant_air_continuum! end

soil_plant_air_continuum!(spac::MonoMLTreeSPAC{FT}, δt::FT; update::Bool = false) where {FT<:AbstractFloat} = (
    # 1. run canopy RT
    canopy_radiation!(spac);

    # 2. run plant hydraulic model
    xylem_pressure_profile!(spac; update = update);

    # 3. run photosynthesis model
    leaf_photosynthesis!(spac, GCO₂Mode());

    # 4. run soil water budget (TODO: add top soil evaporation)
    soil_water!(spac);

    # 5. run soil energy budget
    soil_energy!(spac);

    # 6. run leaf stomatal conductance
    stomatal_conductance!(spac);

    # 7. run leaf energy budget
    # place holder for the function (in SPAC module)

    # 8. update the prognostic variables
    # update soil water content
    # update soil temperature
    stomatal_conductance!(spac, δt);
    # update leaf temperature

    # 9. update xylem flow profiles from stomatal conductance
    xylem_flow_profile!(spac, δt);

    return nothing
);

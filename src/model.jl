#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-12: add function to run soil plant air continuum at a time step
#
#######################################################################################################################################################################################################
"""
This function runs the model using the following steps:
- Run canopy RT module and update the radiation profiles
- Run the photosynthesis and hydraulic models and update the profiles
- Update soil water budget
- Update soil energy budget
- Update leaf stomtal conductances
- Update leaf energy budget

This function is supposed to have the highest hierarchy, and should support all SPAC types defined in ClimaCache.jl

"""
function soil_plant_air_continuum! end

soil_plant_air_continuum!(spac::MonoMLTreeSPAC{FT}) where {FT<:AbstractFloat} = (
    return nothing
);

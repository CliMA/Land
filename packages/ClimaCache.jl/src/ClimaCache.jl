module ClimaCache

using LazyArtifacts

using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using NetcdfIO: read_nc

using EmeraldConstants: CP_D_MOL, CP_L, CP_L_MOL, CP_V_MOL, GAS_R, M_H₂O, P_ATM, T₀, T₂₅, ρ_H₂O, ρg_MPa


#######################################################################################################################################################################################################
#
# Changes to these constants
# General
#     2020-May-30: put the 2017 mat file in an artifact
#     2020-Aug-30: add the updated 2021 mat file along with that of 2017 into a new artifact
#     2022-Jul-20: use reprocessed data for the default constructor
#
#######################################################################################################################################################################################################
const LAND_2017 = artifact"land_model_spectrum_V2" * "/clima_land_spectra_2017.nc";
const LAND_2021 = artifact"land_model_spectrum_V2" * "/clima_land_spectra_2021.nc";


# include the utility types and structures
include("util/colimit.jl")

# include the air types and structures
include("air/air_layer.jl"  )
include("air/meteorology.jl")

# include the radiation types and structures
include("radiation/canopy_optics.jl"      )
include("radiation/canopy_radiation.jl"   )
include("radiation/feature_absorption.jl" )
include("radiation/wave_length_set.jl"    )
include("radiation/canopy.jl"             )
include("radiation/sun_sensor_geometry.jl")
include("radiation/solar_radiation.jl"    )

# include the soil types and structures
include("soil/vulnerability.jl")
include("soil/soil.jl"         )

# include the plant types and structures
include("plant/temperature_dependency.jl")
include("plant/leaf_biophysics.jl"       )
include("plant/leaf_photosynthesis.jl"   )
include("plant/leaf_reaction_center.jl"  )
include("plant/pressure_volume.jl"       )
include("plant/flow.jl"                  )
include("plant/vulnerability.jl"         )
include("plant/hydraulics.jl"            )
include("plant/stomata_model.jl"         )
include("plant/leaf.jl"                  )
include("plant/root.jl"                  )
include("plant/stem.jl"                  )

# include the spac types and structures
include("spac/spac.jl")

end # module

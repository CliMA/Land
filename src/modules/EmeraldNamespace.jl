module EmeraldNamespace

using LazyArtifacts

using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using NetcdfIO: read_nc

using ..EmeraldConstants: CP_D_MOL, CP_L, CP_L_MOL, CP_V_MOL, GAS_R, M_H₂O, P_ATM, T₀, T₂₅, ρ_H₂O, ρg_MPa


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
include("../../packages/EmeraldNamespace.jl/src/util/colimit.jl")

# include the air types and structures
include("../../packages/EmeraldNamespace.jl/src/air/air_layer.jl"  )
include("../../packages/EmeraldNamespace.jl/src/air/meteorology.jl")

# include the radiation types and structures
include("../../packages/EmeraldNamespace.jl/src/radiation/canopy_optics.jl"      )
include("../../packages/EmeraldNamespace.jl/src/radiation/canopy_radiation.jl"   )
include("../../packages/EmeraldNamespace.jl/src/radiation/feature_absorption.jl" )
include("../../packages/EmeraldNamespace.jl/src/radiation/wave_length_set.jl"    )
include("../../packages/EmeraldNamespace.jl/src/radiation/canopy.jl"             )
include("../../packages/EmeraldNamespace.jl/src/radiation/sun_sensor_geometry.jl")
include("../../packages/EmeraldNamespace.jl/src/radiation/solar_radiation.jl"    )

# include the soil types and structures
include("../../packages/EmeraldNamespace.jl/src/soil/vulnerability.jl")
include("../../packages/EmeraldNamespace.jl/src/soil/soil.jl"         )

# include the plant types and structures
include("../../packages/EmeraldNamespace.jl/src/plant/temperature_dependency.jl")
include("../../packages/EmeraldNamespace.jl/src/plant/leaf_biophysics.jl"       )
include("../../packages/EmeraldNamespace.jl/src/plant/leaf_photosynthesis.jl"   )
include("../../packages/EmeraldNamespace.jl/src/plant/leaf_reaction_center.jl"  )
include("../../packages/EmeraldNamespace.jl/src/plant/pressure_volume.jl"       )
include("../../packages/EmeraldNamespace.jl/src/plant/flow.jl"                  )
include("../../packages/EmeraldNamespace.jl/src/plant/vulnerability.jl"         )
include("../../packages/EmeraldNamespace.jl/src/plant/hydraulics.jl"            )
include("../../packages/EmeraldNamespace.jl/src/plant/stomata_model.jl"         )
include("../../packages/EmeraldNamespace.jl/src/plant/leaf.jl"                  )
include("../../packages/EmeraldNamespace.jl/src/plant/root.jl"                  )
include("../../packages/EmeraldNamespace.jl/src/plant/stem.jl"                  )

# include the spac types and structures
include("../../packages/EmeraldNamespace.jl/src/spac/spac.jl")

end # module

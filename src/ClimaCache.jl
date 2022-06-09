#######################################################################################################################################################################################################
#
# Style Guide
# 1. Types and modules start with upper case letters for each word, such as LeafBiophysics
# 2. Functions are lower case words connected using _, such as leaf_biophysics
# 3. Constants are defined using all upper cases, such as LEAF_BIOPHYSICS
# 4. Variables are defined using all lower cases, such as leaf_bio_para
# 5. Temporary variables are defined to start with _, such as _leaf
# 6. Maximum length of each line is 200 letters (including space)
# 7. There should be 2 lines of  empty lines between different components, such as two functions and methods
# 8. Bug fixes or new changes should be documented in the comments above the struct, function, or method, such as this Style Guide above Smaragdus.jl
# 9. Function parameter list that spans multiple lines need to be spaced with 12 spaces (3 tabs)
#
#######################################################################################################################################################################################################
module ClimaCache

using LazyArtifacts

using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using MAT: matread
using PkgUtility: P_ATM, T_0, T_25
using Statistics: mean
using UnPack: @unpack
using WaterPhysics: saturation_vapor_pressure


#######################################################################################################################################################################################################
#
# Changes to these constants
# General
#     2020-May-30: put the 2017 mat file in an artifact
#     2020-Aug-30: add the updated 2021 mat file along with that of 2017 into a new artifact
#
#######################################################################################################################################################################################################
const FILE_SUN    = artifact"land_model_spectrum_V1" * "/sun.mat";
const OPTI_2017   = artifact"land_model_spectrum_V1" * "/Optipar2017_ProspectD.mat";
const OPTI_2021   = artifact"land_model_spectrum_V1" * "/Optipar2021_ProspectPRO_CX.mat";
const WAVELENGTHS = [collect(400:10:650.1); collect(655:5:770.1); collect(780:25:2400.1)];


# include the utility types and structures
include("util/colimit.jl")

# include the air types and structures
include("air/air_layer.jl")

# include the radiation types and structures
include("radiation/canopy_optics.jl"          )
include("radiation/canopy_radiation.jl"       )
include("radiation/canopy.jl"                 )
include("radiation/sun_sensor_geometry.jl"    )
include("radiation/wave_length_set.jl"        )
include("radiation/feature_absorption.jl"     )
include("radiation/hyperspectral_radiation.jl")

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
include("plant/leaf.jl"                  )
include("plant/root.jl"                  )
include("plant/stem.jl"                  )

# include the spac types and structures
include("spac/spac.jl")

end # module

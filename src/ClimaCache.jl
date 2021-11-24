module ClimaCache

using LazyArtifacts

using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using MAT: matread
using Statistics: mean
using UnPack: @unpack


# export public types and structures
export HyperspectralRadiation, WaveLengthSet
export BrooksCorey, VanGenuchten


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


# include the radiation types and structures
include("radiation/wave_length_set.jl")
include("radiation/hyperspectral_radiation.jl")

# include the plant types and structures
include("plant/leaf.jl")

# include the soil types and structures
include("soil/vulnerability_curve.jl")


end # module

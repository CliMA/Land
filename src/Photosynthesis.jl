module Photosynthesis

using ClimaCache: Arrhenius, ArrheniusPeak, Q10
using DocStringExtensions: METHODLIST
using PkgUtility: GAS_R
using UnPack: @unpack


# include the functions
include("temperature.jl")


# include the old code here for back compatibility
# if one need to use the old code, do `using Photosynthesis.PhotosynthesisOld``
include("old/PhotosynthesisOld.jl")


end # module

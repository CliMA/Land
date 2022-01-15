module Photosynthesis

using ClimaCache: AirLayer, Arrhenius, ArrheniusPeak, C3VJPModel, C4VJPModel, GCO₂Mode, Leaf, PCO₂Mode, PhotosynthesisReactionCenter, Q10, VanDerTolFluorescenceModel
using DocStringExtensions: METHODLIST
using PkgUtility: GAS_R, lower_quadratic
using UnPack: @unpack


# export function from this module
export leaf_fluorescence!, leaf_photosynthesis!

# export types from ClimaCache
export AirLayer, GCO₂Mode, Leaf, PCO₂Mode


# include the functions
include("etr.jl")
include("fluorescence.jl")
include("photosynthesis.jl")
include("temperature.jl")

include("old/PhotosynthesisOld.jl")

end # module

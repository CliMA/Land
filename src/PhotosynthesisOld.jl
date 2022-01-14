module PhotosynthesisOld

using ClimaCache: AbstractTemperatureDependency, Arrhenius, ArrheniusPeak, Q10
using DocStringExtensions: TYPEDFIELDS
using PkgUtility: GAS_R, RT_25, T_25, lower_quadratic
using UnPack: @unpack
using WaterPhysics: saturation_vapor_pressure

using ..Photosynthesis: temperature_correction


# export public types
export AirLayer, C3ParaSet, C4ParaSet, GCO₂Mode, Leaf, PCO₂Mode

# export pre-defined parasets
export C3CLM, C4CLM

# export functions
export leaf_fluorescence!, leaf_photosynthesis!, leaf_temperature_dependence!




include("old/types.jl")
include("old/parasets.jl")

include("old/temperature/dependency.jl")

include("old/photosynthesis/etr.jl"           )
include("old/photosynthesis/lightlimited.jl"  )
include("old/photosynthesis/productlimited.jl")
include("old/photosynthesis/rubiscolimited.jl")
include("old/photosynthesis/model.jl"         )

include("old/fluorescence/fluorescence.jl")




end # module

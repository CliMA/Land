module Photosynthesis

using DocStringExtensions: TYPEDFIELDS
using PkgUtility: GAS_R, RT_25, T_25, lower_quadratic
using UnPack: @unpack
using WaterPhysics: saturation_vapor_pressure




# export public types
export AirLayer, C3ParaSet, C4ParaSet, GCO₂Mode, Leaf, PCO₂Mode

# export parasets
export C3Bernacchi, C3CLM, C4CLM

# export functions
export leaf_fluorescence!, leaf_photosynthesis!, leaf_temperature_dependence!




include("types/environment.jl" )
include("types/fluorescence.jl")
include("types/leaf.jl"        )
include("types/mode.jl"        )
include("types/temperature.jl" )
include("types/photomodel.jl"  )
include("types/parasets.jl"    )

include("temperature/correction.jl")
include("temperature/dependency.jl")

include("photosynthesis/etr.jl"           )
include("photosynthesis/lightlimited.jl"  )
include("photosynthesis/productlimited.jl")
include("photosynthesis/rubiscolimited.jl")
include("photosynthesis/model.jl"         )

include("fluorescence/fluorescence.jl")




end # module

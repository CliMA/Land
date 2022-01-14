module PhotosynthesisOld

using DocStringExtensions: TYPEDFIELDS
using PkgUtility: GAS_R, RT_25, T_25, lower_quadratic
using UnPack: @unpack
using WaterPhysics: saturation_vapor_pressure


include("types.jl")
include("parasets.jl")

include("fluorescence/fluorescence.jl")


end # module

module PhotosynthesisOld

using ClimaCache: AbstractTemperatureDependency, Arrhenius, ArrheniusPeak, JmaxTDBernacchi, JmaxTDCLM, KcTDBernacchi, KcTDCLM, KoTDBernacchi, KoTDCLM, KpepTDCLM, Q10, RespirationTDBernacchi,
      RespirationTDCLM, VcmaxTDBernacchi, VcmaxTDCLM, VpmaxTDBoyd, ΓStarTDBernacchi, ΓStarTDCLM
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


include("types.jl")
include("parasets.jl")

include("temperature/dependency.jl")

include("photosynthesis/etr.jl"           )
include("photosynthesis/lightlimited.jl"  )
include("photosynthesis/productlimited.jl")
include("photosynthesis/rubiscolimited.jl")
include("photosynthesis/model.jl"         )

include("fluorescence/fluorescence.jl")


end # module

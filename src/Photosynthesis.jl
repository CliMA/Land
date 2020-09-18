module Photosynthesis

using CLIMAParameters
using CLIMAParameters.Planet
using DocStringExtensions
using Parameters
using WaterPhysics




# define constants here
struct EarthParameterSet <: AbstractEarthParameterSet end
const EARTH = EarthParameterSet()

T_25(FT)  = FT(T_freeze(EARTH)) + 25;
RT_25(FT) = FT(gas_constant()) * T_25(FT);




# export public types
export AbstractFluoModelParaSet,
       AbstractPhotoModelParaSet,
       AbstractTDParameterSet,
       AirLayer,
       ArrheniusPeakTD,
       ArrheniusTD,
       C3ParaSet,
       C4ParaSet,
       FluoParaSet,
       Leaf




# export parasets
export C3Bernacchi,
       C3CLM,
       C4CLM,
       JmaxTDBernacchi,
       JmaxTDCLM,
       JmaxTDLeuning,
       KcTDBernacchi,
       KcTDCLM,
       KoTDBernacchi,
       KoTDCLM,
       KpepTDBoyd,
       KpepTDCLM,
       RespirationTDBernacchi,
       RespirationTDCLM,
       VcmaxTDBernacchi,
       VcmaxTDCLM,
       VcmaxTDLeuning,
       VomaxTDBernacchi,
       VpmaxTDBoyd,
       ΓStarTDBernacchi,
       ΓStarTDCLM




# export functions
export temperature_correction,
       leaf_ETR!,
       leaf_fluorescence!,
       leaf_jmax!,
       leaf_kc!,
       leaf_km!,
       leaf_ko!,
       leaf_kpep!,
       leaf_photo_from_glc!,
       leaf_photo_from_pi!,
       leaf_rd!,
       leaf_temperature_dependence!,
       leaf_vcmax!,
       leaf_vpmax!,
       leaf_Γstar!,
       light_limited_rate!,
       light_limited_rate_glc!,
       photo_TD_from_set,
       photo_TD_from_val,
       product_limited_rate!,
       product_limited_rate_glc!,
       rubisco_limited_rate!,
       rubisco_limited_rate_glc!




include("math/math.jl")

include("types/environment.jl" )
include("types/fluorescence.jl")
include("types/leaf.jl"        )
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

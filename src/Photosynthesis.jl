module Photosynthesis

using CLIMAParameters
using CLIMAParameters.Planet
using DocStringExtensions
using Parameters
using PkgUtility
using WaterPhysics




# define constants here
struct EarthParameterSet <: AbstractEarthParameterSet end
const EARTH = EarthParameterSet()

T_25(FT=Float64)  = FT(T_freeze(EARTH)) + 25;
RT_25(FT=Float64) = FT(gas_constant()) * T_25(FT);




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
       Leaf,
       Q10TD

# export parasets
export C3Bernacchi,
       C3CLM,
       C4CLM,
       FluorescenceVanDerTol,
       FluorescenceVanDerTolDrought,
       JmaxTDBernacchi,
       JmaxTDCLM,
       JmaxTDLeuning,
       KcTDBernacchi,
       KcTDCLM,
       KoTDBernacchi,
       KoTDCLM,
       KpepTDBoyd,
       KpepTDCLM,
       Q10TDAngiosperm,
       Q10TDGymnosperm,
       RespirationTDBernacchi,
       RespirationTDCLM,
       VcmaxTDBernacchi,
       VcmaxTDCLM,
       VcmaxTDLeuning,
       VomaxTDBernacchi,
       VpmaxTDBoyd,
       VtoRCollatz,
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
       leaf_photosynthesis!,
       leaf_rd!,
       leaf_temperature_dependence!,
       leaf_vcmax!,
       leaf_vpmax!,
       leaf_Γstar!,
       light_limited_rate!,
       photo_TD_from_set,
       photo_TD_from_val,
       product_limited_rate!,
       rubisco_limited_rate!




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

module Photosynthesis

using CLIMAParameters
using DocStringExtensions
using Parameters
using WaterPhysics

# define constants here
const GAS_R = gas_constant()
const K_25  = 298.15

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

# export functions
export arrhenius_correction,
       leaf_ETR!,
       leaf_ETR_Jps,
       leaf_ETR_pot_APAR,
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
       light_limited_an_glc,
       light_limited_rate!,
       light_limited_rate_glc!,
       photo_TD_from_set,
       photo_TD_from_val,
       product_limited_an_glc,
       product_limited_rate!,
       product_limited_rate_glc!,
       rubisco_limited_an_glc,
       rubisco_limited_rate!,
       rubisco_limited_rate_glc!




include("types.jl"       )
include("parasets.jl"    )
include("math.jl"        )
include("temperature.jl" )
include("photorates.jl"  )
include("photomodel.jl"  )
include("fluorescence.jl")




end # module

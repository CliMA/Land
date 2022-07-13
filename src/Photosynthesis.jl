module Photosynthesis

using ClimaCache: AirLayer, AndereggSM, Arrhenius, ArrheniusPeak, BallBerrySM, BetaFunction, BetaParameterG1, BetaParameterVcmax, C3VJPModel, C3CytochromeModel, C4VJPModel, CytochromeReactionCenter,
      EllerSM, GentineSM, GCO₂Mode, Leaf, Leaves1D, Leaves2D, LeuningSM, MedlynSM, MinimumColimit, MonoElementSPAC, MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC, PCO₂Mode, Q10, QuadraticColimit,
      SerialColimit, SperrySM, VJPReactionCenter, VanDerTolFluorescenceModel, WangSM, Wang2SM
using PkgUtility: GAS_R, lower_quadratic, upper_quadratic
using UnPack: @unpack


# export function from this module
export leaf_photosynthesis!


# include the functions
include("colimit.jl"        )
include("etr.jl"            )
include("fluorescence.jl"   )
include("light_limited.jl"  )
include("model.jl"          )
include("product_limited.jl")
include("rubisco_limited.jl")
include("temperature.jl"    )


end # module

module Photosynthesis

using ClimaCache: CytochromeReactionCenter, VJPReactionCenter, VanDerTolFluorescenceModel
using ClimaCache: Arrhenius, ArrheniusPeak, C3VJPModel, C3CytochromeModel, C4VJPModel, GCO₂Mode, MinimumColimit, PCO₂Mode, Q10, QuadraticColimit, SerialColimit
using ClimaCache: AndereggSM, BallBerrySM, BetaFunction, BetaParameterG1, BetaParameterVcmax, EllerSM, GentineSM, LeuningSM, MedlynSM, SperrySM, WangSM, Wang2SM
using ClimaCache: AirLayer, Leaf, Leaves1D, Leaves2D
using ClimaCache: MonoElementSPAC, MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC
using ClimaCache: F_O₂, GAS_R
using PkgUtility: lower_quadratic, upper_quadratic
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

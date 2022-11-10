module Photosynthesis

using PkgUtility: lower_quadratic, upper_quadratic
using UnPack: @unpack

using EmeraldConstants: F_O₂, GAS_R
using ClimaCache: CytochromeReactionCenter, VJPReactionCenter, VanDerTolFluorescenceModel
using ClimaCache: Arrhenius, ArrheniusPeak, C3VJPModel, C3CytochromeModel, C4VJPModel, GCO₂Mode, PCO₂Mode, Q10, Q10Peak
using ClimaCache: MinimumColimit, QuadraticColimit, SerialColimit, SquareColimit
using ClimaCache: AndereggSM, BallBerrySM, BetaFunction, BetaParameterG1, BetaParameterVcmax, EllerSM, GentineSM, LeuningSM, MedlynSM, SperrySM, WangSM, Wang2SM
using ClimaCache: AirLayer, Leaf, Leaves1D, Leaves2D
using ClimaCache: MonoElementSPAC, MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC


include("colimit.jl"        )
include("etr.jl"            )
include("fluorescence.jl"   )
include("light_limited.jl"  )
include("model.jl"          )
include("product_limited.jl")
include("rubisco_limited.jl")
include("temperature.jl"    )


end # module

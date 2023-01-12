module Photosynthesis

using PkgUtility: lower_quadratic, upper_quadratic
using UnPack: @unpack

#using EmeraldConstants: F_O₂, GAS_R
#using EmeraldNamespace: CytochromeReactionCenter, VJPReactionCenter, VanDerTolFluorescenceModel
#using EmeraldNamespace: Arrhenius, ArrheniusPeak, C3VJPModel, C3CytochromeModel, C4VJPModel, GCO₂Mode, PCO₂Mode, Q10, Q10Peak
#using EmeraldNamespace: MinimumColimit, QuadraticColimit, SerialColimit, SquareColimit
#using EmeraldNamespace: AbstractStomataModel, BallBerrySM, BetaFunction, BetaParameterG1, BetaParameterVcmax, GentineSM, LeuningSM, MedlynSM
#using EmeraldNamespace: AirLayer, Leaf, Leaves1D, Leaves2D
#using EmeraldNamespace: MonoElementSPAC, MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC


include("colimit.jl"        )
include("etr.jl"            )
include("fluorescence.jl"   )
include("light_limited.jl"  )
include("model.jl"          )
include("product_limited.jl")
include("rubisco_limited.jl")
include("temperature.jl"    )


end # module

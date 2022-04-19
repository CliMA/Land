module SoilHydraulics

using ClimaCache: VanGenuchten
using ConstrainedRootSolvers: ReduceStepMethodND, SolutionToleranceND, find_peak
using DocStringExtensions: METHODLIST
using UnPack: @unpack

import ClimaCache: BrooksCorey


include("vulnerability.jl")


end # module

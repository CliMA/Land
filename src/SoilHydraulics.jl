module SoilHydraulics

using ClimaCache: VanGenuchten
using ConstrainedRootSolvers: ReduceStepMethodND, SolutionToleranceND, find_peak
using DocStringExtensions: METHODLIST
using UnPack: @unpack

import ClimaCache: BrooksCorey


export BrooksCorey, VanGenuchten, relative_hydraulic_conductance,  soil_ψ_25, soil_θ


include("vulnerability.jl")


end # module

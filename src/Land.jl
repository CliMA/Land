module Land

using DiffEqOperators
using DocStringExtensions
using OrdinaryDiffEq
using Parameters
using WaterPhysics: TraceGasAir, TraceGasCO₂, diffusive_coefficient


# include types
include("Land/Types.jl"    )
include("Land/AirLayers.jl")


end # module

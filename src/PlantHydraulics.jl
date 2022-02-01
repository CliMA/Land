module PlantHydraulics

using ClimaCache: ComplexVC, LogisticVC, PowerVC, WeibullVC
using DocStringExtensions: METHODLIST
using UnPack: @unpack


# export public types from ClimaCache
export ComplexVC, LogisticVC, PowerVC, WeibullVC

# export public functions
export relative_hydraulic_conductance


# include functions
include("vulnerability.jl")

# include old module
include("old/PlantHydraulicsOld.jl")


end # module

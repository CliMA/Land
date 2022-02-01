module PlantHydraulics

using ClimaCache: LogisticVC, PowerVC, WeibullVC
using DocStringExtensions: METHODLIST
using UnPack: @unpack


# export public types from ClimaCache
export LogisticVC, PowerVC, WeibullVC

# export public functions
export xylem_k_ratio


# include functions
include("vulnerability.jl")

# include old module
include("old/PlantHydraulicsOld.jl")


end # module

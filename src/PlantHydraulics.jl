module PlantHydraulics

using ClimaCache: ComplexVC, LogisticVC, PowerVC, WeibullVC
using DocStringExtensions: METHODLIST
using UnPack: @unpack

import SoilHydraulics: relative_hydraulic_conductance


# export public types from ClimaCache
export ComplexVC, LogisticVC, PowerVC, WeibullVC

# export public functions
export critical_pressure, relative_hydraulic_conductance


# include functions
include("critical_pressure.jl")
include("vulnerability.jl"    )

# include old module
include("old/PlantHydraulicsOld.jl")


end # module

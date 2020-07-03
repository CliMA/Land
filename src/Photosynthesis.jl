module Photosynthesis

using CLIMAParameters
using DocStringExtensions
using Parameters

# define constants here
const GAS_R = gas_constant()
const K_25  = 298.15




include("types.jl"     )
include("parasets.jl"  )
include("math.jl"      )
include("photomodel.jl")




end # module

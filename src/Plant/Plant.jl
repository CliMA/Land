module Plant

using DocStringExtensions
using Parameters

using ..LandParameters
using ..WaterPhysics
using ..Photosynthesis

@unpack GRAVITY,
        K_25,
        ρ_H₂O = LandParameters

export AbstractEmpiricalStomatalModel,
       AbstractOptimizationStomatalModel,
       ESMBallBerry,
       Tree,
       get_empirical_gsw_pi




include("types.jl")
include("hydraulics.jl")
include("stomata.jl")




###############################################################################
#
# Functions in development
# Test and document the functions before merging them to this file
#
###############################################################################
include("Plant_in_development.jl")




end

module Plant

#=
using DocStringExtensions
using Parameters
using RootSolvers

using ..LandParameters
using ..WaterPhysics
using ..Photosynthesis

@unpack GRAVITY,
        K_25,
        ρ_H₂O = LandParameters




include("types.jl"     )
include("hydraulics.jl")
include("stomata.jl"   )
include("testing.jl"   )




###############################################################################
#
# Functions in development
# Test and document the functions before merging them to this file
#
###############################################################################
include("Plant_in_development.jl")
=#



end

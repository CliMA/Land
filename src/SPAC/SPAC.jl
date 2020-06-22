module SPAC

#=
using DocStringExtensions
using Parameters

using ..WaterPhysics
using ..Plant
using ..CanopyRT




include("types.jl"      )
include("paraset.jl"    )
include("simulations.jl")




###############################################################################
#
# Functions in development
# Test and document the functions before merging them to this file
#
###############################################################################
include("SPAC_in_development.jl")
=#



end # module

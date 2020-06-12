module Leaf

using DocStringExtensions
using Parameters

using ..LandParameters
using ..WaterPhysics
using ..Photosynthesis
using ..Plant

@unpack CP_D,
        GRAVITY,
        R_D,
        VON_KARMAN_CONST,
        WATER_AIR_MRATIO = LandParameters

export LeafParams,
       MeteoParams,
       leaf_photosynthesis!




###############################################################################
#
# Types in development
# Test and document the functions before merging them to this file
#
###############################################################################
include("Type_in_development.jl")




include("types.jl"    )
include("parasets.jl" )
include("leafphoto.jl")




###############################################################################
#
# Functions in development
# Test and document the functions before merging them to this file
#
###############################################################################
include("Leaf_in_development.jl")




end
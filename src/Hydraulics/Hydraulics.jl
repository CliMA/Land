module Hydraulics

using DocStringExtensions
using Parameters
using RootSolvers

using ..MathTools

# TODO add extra-xylary VC
# TODO add pressure-volume curve functions

# export public types
export LeafHydraulics,
       StemHydraulics

# export public functions
export hydraulic_p_profile!,
       leaf_e_crit,
       leaf_xylem_risk,
       xylem_p_from_flow




include("types.jl")
include("flow.jl" )




end # module

module Hydraulics

using DocStringExtensions
using Parameters
using RootSolvers

using ..LandParameters
using ..MathTools
using ..RootSolversExtension

# TODO add extra-xylary VC
# TODO add pressure-volume curve functions
# TODO Bridge to a more general vG formulations

const GRAVITY = LandParameters.GRAVITY
const K_25    = LandParameters.K_25
const ρ_H₂O   = LandParameters.ρ_H₂O
const ρg_MPa  = ρ_H₂O * GRAVITY * 1e-6

# export public types
export LeafHydraulics,
       RootHydraulics,
       StemHydraulics,
       Tree

# export public functions
export create_tree,
       hydraulic_p_profile!,
       leaf_e_crit,
       leaf_xylem_risk,
       root_qs_p_from_q,
       xylem_p_from_flow




include("types.jl"  )
include("testing.jl")
include("flow.jl"   )




end # module
